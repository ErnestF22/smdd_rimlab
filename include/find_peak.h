#ifndef FIND_PEAK_H
#define FIND_PEAK_H

#include <rofl/common/macros.h>
#include <vector>

inline unsigned int index_diff(unsigned int i1,
                               unsigned int i2,
                               unsigned int n)
{
    if (i1 < i2)
        return ((i2 - i1) % n);
    else
        return ((i1 - i2) % n);
}

template <typename Func, typename T, typename InserIt>
void find_peak(Func func,
               unsigned int ibeg,
               unsigned int iend,
               unsigned int win,
               const T &minval,
               InserIt ins)
{
    unsigned i, imax;
    unsigned j, jbeg, jend, jmax;
    std::vector<unsigned int> candidates;
    // Visits the indices from ibeg to iend and find maxima of function func
    // (which maps an unsigned integer to a T value).
    // A maximum M is such that func(i) <= func(M) for M-win <= i <= M+win.
    i = ibeg;
    imax = ibeg;
    while (i < iend)
    {
        // std::cout << "item i: " << i << " f(i) " << func(i) << ", imax " << imax
        //   << " f(imax) " << func(imax) << "\n";
        jmax = i;
        // If the previous maximum is not in the interval, then the new imax
        // must be computed on the whole window [i-win,i+win] and not only on
        // [i,i+win]
        if (imax + win < i)
        {
            jbeg = (i >= win ? i - win : 0);
            imax = i;
        }
        else
        {
            jbeg = i;
        }
        jend = (i + win + 1 < iend ? i + win + 1 : iend);
        // std::cout << "  visiting items j:" << jbeg << "," << jend << "\n";
        for (j = jbeg; j < jend; ++j)
        {
            if (func(j) > func(jmax))
            {
                jmax = j;
            }
        }
        // The new imax is the maximum between the previous imax (if valid...)
        // and jmax
        imax = (func(jmax) >= func(imax) ? jmax : imax);
        // std::cout << "  jmax " << jmax << " f(jmax) " << func(jmax)
        //           << ", new imax " << imax << "\n ";
        if (i == imax && func(i) > minval)
        {
            // std::cout << "  *** local peak in " << i << ", f(i) " << func(i)
            //           << "\n";
            candidates.push_back(i);
        }
        // Advance i
        if (jmax > i)
            i = jmax;
        else
            ++i;
    }

    // Removes peaks that lie at a distance less than win
    for (int i = 0; i < (int)candidates.size() - 1; ++i)
    {
        ROFL_VAR2(i, (int)candidates.size())
        if (candidates[i] - candidates[i + 1] < win)
            if (candidates[i] >= candidates[i + 1])
            {
                candidates.erase(candidates.begin() + i + 1);
                i--;
            }
            else
            {
                candidates.erase(candidates.begin() + i);
                i--;
            }
    }

    for (auto &c : candidates)
        ins = c;

    // std::vector<unsigned int>::iterator it, mid;
    // int counter = 0;
    // it = candidates.begin();
    // mid = it;
    // for (it = candidates.begin(); it != candidates.end(); ++it) {
    //     ROFL_VAR4(*it, *mid, func(*mid), std::distance(mid, it));
    //     //    std::cerr << "*it " << *it << " val " << func(*it)
    //     //      << ", *mid " << *mid << " val " << func(*mid)
    //     //      << ", distance(mid,it) " << std::distance(mid,it) << "\n";
    //     if (func(*mid) == func(*it)) {
    //         counter++;
    //         if (2 * std::distance(mid, it) > counter) {
    //             ++mid;
    //         }
    //         std::cout << "  counter " << counter << ", *mid " << *mid
    //                   << std::endl;
    //     } else {
    //         ins = *mid;
    //         std::cout << "  save peak " << *mid << " value " << func(*mid)
    //                   << std::endl;
    //         mid = it;
    //         counter = 0;
    //     }
    // }
    // if (mid != candidates.end()) {
    //     ins = *mid;
    // }
}

template <typename Func, typename T, typename InserIt>
void find_peak_circular(Func func,
                        unsigned int ibeg,
                        unsigned int iend,
                        unsigned int win,
                        const T &minval,
                        InserIt ins)
{
    // Finds all maxima
    std::vector<int> maxima;
    find_peak(func, ibeg, iend, win, minval, std::back_inserter(maxima));

    // Keeps only the maxima that do not lies in a window between the beginning
    // and end of interval
    unsigned int nmax = maxima.size();
    unsigned int period = iend - ibeg;
    bool keepFirst = true;
    bool keepLast = true;
    // if (nmax > 1 && index_diff(maxima[nmax - 1], maxima[0], nmax) <= win) {
    //  if (nmax > 1 && (maxima[0] + period - maxima[nmax - 1]) <= win) {
    //      if (func(maxima[0]) < func(maxima[nmax - 1])) {
    //          keepFirst = false;
    //      } else if (func(maxima[0]) > func(maxima[nmax - 1])) {
    //          keepLast = false;
    //      }
    //  }
    if (maxima.empty())
        return;

    if (maxima[0] < win)
        for (unsigned int i = 1; i <= win && keepFirst; ++i)
        {
            if (func(maxima[0]) < func((maxima[0] + period - i) % period))
            {
                keepFirst = false;
                ROFL_VAR1("keepFirst false")
                break;
            }
        }
    if (maxima.back() > period - win)
        for (unsigned int i = 1; i <= win && keepLast; ++i)
        {
            if (func(maxima.back()) < func((maxima.back() + i) % period))
            {
                // some comparisons are redundant
                keepLast = false;
                ROFL_VAR6(maxima.back(), (maxima.back() + i) % period, func(maxima.back()), func((maxima.back() + i) % period), i, "keepLast false")
                break;
            }
        }

    // Inserts the value of maxima that are not shadowed by another maximum
    for (unsigned int i = 0; i < nmax; ++i)
    {
        if ((i != 0 && i != nmax - 1) || (i == 0 && keepFirst) ||
            (i == nmax - 1 && keepLast))
        {
            ins = maxima[i];
        }
    }
}

// template <typename Cmp,typename InserIt>
// void find_peak2(unsigned int ibeg,unsigned int iend,unsigned int win,
//                  Cmp cmp,bool circular,InserIt ins)
//{
//   assert(ibeg <= iend);
//   unsigned int n = iend - iend;
//   std::vector<unsigned int> heap;  // to keep the maximum value

//  std::cout << "ibeg " << ibeg << ", iend " << iend << ", win " << win <<
//  std::endl;

//  // Initializes heap with item inside a window around ibeg
//  heap.push_back(ibeg);
//  for (unsigned int j = 1; j < n && j <= win; ++j) {
//    heap.push_back(ibeg + j);
//    if (circular && win < n) {
//      heap.push_back(ibeg + n - j);
//    }
//  }
//  std::make_heap(heap.begin(),heap.end(),cmp);

//  std::cout << "heap: ";
//  for (int h = 0; h < heap.size(); ++h)
//    std::cout << heap[h] << " ";
//  std::cout << std::endl;
//
//  // Visits all the histogram and pushes the items into a heap
//  for (unsigned int i = ibeg; i < iend; ++i) {
//    // Extracts values from the heap while they do not lie in the window
//    around i while (!heap.empty() && index_diff(heap.front(),i,n) > win) {
//      std::pop_heap (heap.begin(),heap.end());
//      heap.pop_back();
//    }
//    // Check if i is the maximum inside the heap
//    if (!heap.empty() && cmp(i,heap.front()) == true) {
//      ins = i;
//    }
//    // Inserts the value of the new end point of window inside the heap
//    if (!circular && i+win+1 < iend) {
//      heap.push_back(i + 1 + win);
//      std::push_heap(heap.begin(),heap.end(),cmp);
//    }
//    else if (circular) {
//      unsigned int v = (i + 1 + win - ibeg) % n + ibeg;
//      heap.push_back(v);
//      std::push_heap(heap.begin(),heap.end(),cmp);
//    }
//  }
//}

// template <typename Func,typename T,typename InserIt>
// void find_peak2(Func func,unsigned int ibeg,unsigned int iend,unsigned int
// win,
//                 const T& minval,InserIt ins)
//{
//   // It requires that ibeg > iend, otherwisw it returns an empty set
//   if (ibeg > iend) return;

//  // Two vectors defined as follows
//  // - fwd[i]: contains the index of maximum index
//  unsigned int n = iend - ibeg;
//  std::vector<unsigned int> fwd(n,minVal);
//  std::vector<unsigned int> bwd(n,minVal);

//  for (unsigned int i = 0; i < n; ++i) {
//    if (i % win == 0) {
//      fwd[i] = i+ibeg;
//    }
//    else {
//      if (func(i+ibeg) > fwd[i-1]) fwd[i] = i + ibeg;
//      else fwd[i] = fwd[i-1];
//    }
//  }
//  for (unsigned int i = n-1; i >= 0; --i) {
//    if (i == n-1 || i % win == 0) {
//      bwd[i] = i+ibeg;
//    }
//    else {
//      if (func(i+ibeg) > bwd[i+1]) bwd[i] = i + ibeg;
//      else bwd[i] = bwd[i+1];
//    }
//  }
//}

#endif
