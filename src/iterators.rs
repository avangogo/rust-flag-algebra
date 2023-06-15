//! Streaming iterators on combinatorial objects (subsets, functions, ...).
//!
//! All these iterators works without allocating new memory after
//! their initialisation.
///
/// Interface for streaming iterators.
///
/// Similarly as in the streaming-iterator crate, the elements yielded by the iterator
/// are borrowed by the iterator.
/// A loop on such an iterator `iter` is written as follows.
///```ignore
///while let Some(item) = iter.next() {
///    ...
///}
///```

pub trait StreamingIterator<A>
where
    A: ?Sized,
{
    /// Return the next value of the iterator.
    fn next(&mut self) -> Option<&A>;

    /// Consume the iterator and return the number of elements yielded.
    fn count(mut self) -> usize
    where
        Self: Sized,
    {
        let mut count = 0;
        while self.next().is_some() {
            count += 1;
        }
        count
    }
}

/// Iterator on subsets of `[n]`.
#[derive(Clone, Debug)]
pub struct Subsets {
    n: usize,
    data: Vec<usize>,
    untouched: bool,
}

impl Subsets {
    /// Create an iterator on the subsets of `[n]`.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            data: Vec::with_capacity(n),
            untouched: true,
        }
    }
}

impl StreamingIterator<[usize]> for Subsets {
    fn next(&mut self) -> Option<&[usize]> {
        if self.untouched {
            //skip the first step
            self.untouched = false;
        } else {
            let mut k = self.n;
            loop {
                if k == 0 {
                    return None;
                }
                k -= 1;
                if self.data.last() != Some(&k) {
                    break;
                }
                let _ = self.data.pop();
            }
            self.data.push(k);
        }
        Some(&self.data)
    }
}

/// Iterator on functions from [n] to [k].
#[derive(Clone, Debug)]
pub struct Functions {
    n: usize,
    k: usize,
    data: Vec<usize>,
    untouched: bool,
}

impl Functions {
    /// Create an iterator on the function from `[n]` to `[k]`.
    pub fn new(n: usize, k: usize) -> Self {
        Self {
            n,
            k,
            data: vec![0; n],
            untouched: true,
        }
    }
}

impl StreamingIterator<[usize]> for Functions {
    fn next(&mut self) -> Option<&[usize]> {
        if self.untouched {
            //skip the first step
            self.untouched = false;
        } else {
            let mut i = 0;
            while i < self.n && self.data[i] == self.k - 1 {
                i += 1;
            }
            if i == self.n {
                return None;
            } else {
                self.data[i] += 1;
                for j in 0..i {
                    self.data[j] = 0;
                }
            }
        }
        Some(&self.data)
    }
}

/// Iterator on subsets of [n] with `k` elements represented by
/// an increasing array.
#[derive(Clone, Debug)]
pub struct Choose {
    n: usize,     // total number of elements
    k: usize,     // number of elements chosen
    fixed: usize, // number of elements fixed
    data: Vec<usize>,
    untouched: bool,
}

impl Choose {
    /// Subsets of [n] with k elements.
    pub fn new(n: usize, k: usize) -> Self {
        Self::with_fixed_part(n, k, 0)
    }
    /// Subsets of [n] with k elements that contain [fixed].  
    pub fn with_fixed_part(n: usize, k: usize, fixed: usize) -> Self {
        assert!(fixed <= k && k <= n);
        Self {
            n,
            k,
            fixed,
            data: (0..k).collect(),
            untouched: true,
        }
    }
}

impl StreamingIterator<[usize]> for Choose {
    fn next(&mut self) -> Option<&[usize]> {
        if self.untouched {
            self.untouched = false;
        } else {
            let mut i = self.k;
            loop {
                if i <= self.fixed {
                    return None;
                } else {
                    i -= 1
                };
                if self.data[i] != self.n - self.k + i {
                    break;
                }
            }
            self.data[i] += 1;
            for j in (i + 1)..self.k {
                self.data[j] = self.data[i] + j - i;
            }
        }
        Some(&self.data)
    }
}

/// Iterator on the partitions of [n] into two sets of respective size `k` and `n-k`.
#[derive(Clone, Debug)]
pub struct Split {
    choose: Choose,
    second_part: Vec<usize>,
}

impl Split {
    /// Create an iterator on partitions of [n] into two sets of respective size `k` and `n-k`.
    #[allow(unused)]
    pub fn new(n: usize, k: usize) -> Self {
        Self::with_fixed_part(n, k, 0)
    }
    /// Create an iterator on partitions of [n] into two sets of respective size `k`
    /// and `n-k+fixed` that both contain `[fixed]`.
    pub fn with_fixed_part(n: usize, k: usize, fixed: usize) -> Self {
        let second_part_size = n - k + fixed;
        Self {
            choose: Choose::with_fixed_part(n, k, fixed),
            second_part: (0..second_part_size).collect(),
        }
    }
    /// Yield next element of the iterator.
    ///
    /// Because of the type of this function, `Split` does not implement
    /// the trait `StreamingIterator`.
    pub fn next(&mut self) -> Option<(&[usize], &[usize])> {
        let fixed = self.choose.fixed;
        let k = self.choose.k;
        let n = self.choose.n;
        match self.choose.next() {
            Some(p) => {
                let mut j = fixed;
                for i in fixed..=k {
                    let start = if i == fixed { fixed } else { p[i - 1] + 1 };
                    let end = if i == k { n } else { p[i] };
                    for v in start..end {
                        self.second_part[j] = v;
                        j += 1;
                    }
                }
                Some((p, &self.second_part))
            }
            None => None,
        }
    }
}

/// Iterator on the injections from [k] to [n].
#[derive(Clone, Debug)]
pub struct Injection {
    n: usize,
    k: usize,
    data: Vec<usize>,
    index: Vec<usize>,
    fixed: usize, // number of fixed elements
    untouched: bool,
}

impl Injection {
    /// Iterator on the injections from `[k]` to `[n]` that stabilize `[fixed]`.
    pub fn with_fixed_part(n: usize, k: usize, fixed: usize) -> Self {
        assert!(k <= n);
        assert!(fixed <= k);
        Self {
            n,
            k,
            data: (0..n).collect(),
            index: (0..k).collect(),
            fixed,
            untouched: true,
        }
    }

    /// Iterator on the injections from [k] to [n].
    pub fn new(n: usize, k: usize) -> Self {
        Self::with_fixed_part(n, k, 0)
    }

    ///Iterator on the permutations of [n].
    pub fn permutation(n: usize) -> Self {
        Self::new(n, n)
    }

    #[inline]
    fn set_index(&mut self, i: usize, val: usize) {
        let old_val = self.index[i];
        self.index[i] = val;
        self.data.swap(i, old_val);
        self.data.swap(i, val);
    }
}

impl StreamingIterator<[usize]> for Injection {
    fn next(&mut self) -> Option<&[usize]> {
        if self.untouched {
            self.untouched = false;
        } else {
            if self.k == self.fixed {
                return None;
            }; // Avoiding negative value on next line
            let mut i = self.k - 1;
            while self.index[i] == self.n - 1 {
                if i == self.fixed {
                    return None;
                };
                i -= 1
            }
            let vi = self.index[i] + 1;
            self.set_index(i, vi);
            for j in (i + 1)..self.k {
                self.set_index(j, j);
            }
        }
        Some(&self.data[0..self.k])
    }
}

///Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_subsets() {
        assert_eq!(1, Subsets::new(0).count());
        assert_eq!(2, Subsets::new(1).count());
        assert_eq!(16, Subsets::new(4).count());
    }

    #[test]
    fn unit_choose() {
        assert_eq!(120, Choose::new(10, 3).count());
        assert_eq!(3, Choose::new(3, 1).count());
        assert_eq!(1, Choose::new(0, 0).count());
        assert_eq!(1, Choose::new(2, 2).count());
        assert_eq!(120, Choose::with_fixed_part(11, 4, 1).count());
        assert_eq!(1, Choose::with_fixed_part(5, 5, 4).count());
        assert_eq!(1, Choose::with_fixed_part(4, 4, 4).count());
        assert_eq!(1, Choose::with_fixed_part(6, 2, 2).count());
        assert_eq!(1, Choose::with_fixed_part(0, 0, 0).count());
    }

    #[test]
    fn unit_split() {
        let mut iter = Split::with_fixed_part(12, 5, 2);
        let mut count = 0;
        while let Some((_a, _b)) = iter.next() {
            count += 1;
        }
        assert_eq!(120, count);
    }

    #[test]
    fn unit_injection() {
        assert_eq!(60, Injection::new(5, 3).count());
        assert_eq!(1, Injection::new(42, 0).count());
        assert_eq!(720, Injection::permutation(6).count());
        assert_eq!(20, Injection::with_fixed_part(8, 5, 3).count());
        assert_eq!(1, Injection::with_fixed_part(6, 2, 2).count());
        assert_eq!(1, Injection::new(0, 0).count());
    }
}
