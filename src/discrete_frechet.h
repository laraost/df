#ifndef DISCRETE_FRECHET_H
#define DISCRETE_FRECHET_H

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

namespace df {

template<typename distance_func, typename ForwardIterator>
auto compute_discrete_frechet(ForwardIterator p_begin, ForwardIterator p_end,
                              ForwardIterator q_begin, ForwardIterator q_end,
                              distance_func dist_func = {})
{
    using distance_t = decltype(distance_func()(*p_begin, *q_begin));
    using row_t = std::vector<distance_t>;
    using size_t = row_t::size_type;
    std::array<row_t, 2> rows;
    row_t* current_row = &rows[0];
    row_t* next_row = &rows[1];
    const size_t p_length = std::distance(p_begin, p_end);
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
