/*
 * Copyright (c) 2023 Lara Ost
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef DISCRETE_FRECHET_H
#define DISCRETE_FRECHET_H

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
#include <vector>

namespace df {

namespace internal {

//
// Computes the discrete Fréchet distance between a trajectory consisting of a single point
// (given by `p_begin`) and a trajectory Q given by the iterator pair `q_begin`, `q_end`.
//
// Time: O(|Q|)
// Space: O(1)
//
template<typename distance_func, typename ForwardIterator>
auto compute_discrete_frechet_single(ForwardIterator p_point,
                                     ForwardIterator q_begin, ForwardIterator q_end,
                                     distance_func dist_func = {}) {
    auto q_iter = q_begin;
    auto result = dist_func(*p_point, *q_iter);
    for (++q_iter; q_iter != q_end; ++q_iter) {
        result = std::max(result, dist_func(*p_point, *q_iter));
    }
    return result;
}

} // End of namespace `internal`

//
// Computes the discrete Fréchet distance between trajectories P and Q
// using the dynamic programming algorithm by Eiter and Mannila (1994).
//
// Trajectories P and Q are given as pairs of forward iterators:
// P as p_begin and p_end, Q as q_begin and q_end,
// where *_end points one past the last point of the trajectory.
//
// The distance function is given as `dist_func`,
// a binary function of two `decltype(*p_begin)`s and returning a numeric type.
//
// Time: O(|P||Q|)
// Space: O(|P|)
//
template<typename distance_func, typename ForwardIterator>
auto compute_discrete_frechet(ForwardIterator p_begin, ForwardIterator p_end,
                              ForwardIterator q_begin, ForwardIterator q_end,
                              distance_func dist_func = {})
{
    const size_t p_length = std::distance(p_begin, p_end);
    const size_t q_length = std::distance(q_begin, q_end);
    if (p_length == 1) { return internal::compute_discrete_frechet_single(p_begin, q_begin, q_end, dist_func); }
    if (q_length == 1) { return internal::compute_discrete_frechet_single(q_begin, p_begin, p_end, dist_func); }

    using distance_t = decltype(distance_func()(*p_begin, *q_begin));
    using row_t = std::vector<distance_t>;
    using size_t = row_t::size_type;
    std::array<row_t, 2> rows;
    row_t* current_row = &rows[0];
    row_t* next_row = &rows[1];
    current_row->resize(p_length);
    next_row->resize(p_length);

    auto q_iter = q_begin;
    auto p_iter = p_begin;
    // Initialize the first row
    (*current_row)[0] = dist_func(*p_iter, *q_iter);
    ++p_iter;
    for (size_t i = 1; i < p_length; ++i, ++p_iter) {
        const auto d_pq = dist_func(*p_iter, *q_iter);
        (*current_row)[i] = std::max((*current_row)[i-1], d_pq);
    }

    // Compute remaining rows
    for (++q_iter; q_iter != q_end; ++q_iter) {
        p_iter = p_begin;
        (*next_row)[0] = std::max((*current_row)[0], dist_func(*p_iter, *q_iter));
        ++p_iter;
        for (size_t i = 1; i < p_length; ++i, ++p_iter) {
            const auto d_pq = dist_func(*p_iter, *q_iter);
            (*next_row)[i] = std::max(std::min({(*next_row)[i-1],
                                                (*current_row)[i-1],
                                                (*current_row)[i]}),
                                      d_pq);
        }
        std::swap(current_row, next_row);
    }
    return current_row->back();
}

} // End of namespace `df`

#endif
