#ifndef LEMON_PARALLEL_HPP
#define LEMON_PARALLEL_HPP

#include <list>

#include "lemon/hadoop.hpp"
#include "lemon/entries.hpp"

#include <chrono>

#ifdef LEMON_USE_ASYNC
#include "lemon/thread_pool.hpp"
#else
#include <thread>
#endif

#ifdef LEMON_BENCHMARK
#include <iostream>
#endif

namespace lemon {

#ifndef LEMON_USE_ASYNC

//! The `run_parallel` function launches jobs which do return data.
//!
//! Use this function to run the `worker` function on `ncpu` threads. The
//! `worker` should accept two arguments, a `chemfiles::Frame` and a
//! `std::string`. It must return a value as this value will be appended, using
//! the `combine` function object, the the `collector`. See the `Lemon Workflow`
//! documention for more details. \param worker A function object (C++11 lambda,
//! struct the with operator()
//!  overloaded, or std::function object) that the user wishes to apply.
//! \param [in] p A path to the Hadoop sequence file directory.
//!  by the worker object which are appended with `combine`.
//! \param collector A function object that handles the output of `worker`.
//! \param [in] ncpu The number of threads to use.
//! \param [in] entries Which entries to use. Not used if blank.
//! \param [in] skip_entries Which entries to skip. Not used if blank.
template <typename Function, typename Collector>
inline void run_parallel(Function&& worker, const std::string& p,
                         Collector& collector, size_t ncpu = 1,
                         const Entries& entries = Entries(),
                         const Entries& skip_entries = Entries()) {
    auto pathvec = read_hadoop_dir(p);
    std::vector<std::thread> threads(ncpu);
    using ret = typename std::result_of<Function&(chemfiles::Frame,
                                                  const std::string&)>::type;

    // Total number of jobs for each thread
    const int grainsize = static_cast<int>(pathvec.size() / ncpu);
    auto work_iter = pathvec.begin();
    using iter = std::vector<std::string>::iterator;
    std::unordered_map<std::thread::id, std::list<ret>> results;

    auto call_function = [&worker, &results, &entries,
                          &skip_entries](iter first, iter last) {
        auto th = std::this_thread::get_id();
        for (auto it = first; it != last; ++it) {
            std::ifstream data(*it, std::istream::binary);
            Hadoop sequence(data);

            while (sequence.has_next()) {
                auto pair = sequence.next();
                if (entries.size() && entries.count(pair.first) == 0) {
                    continue;
                }
                if (skip_entries.size() &&
                    skip_entries.count(pair.first) != 0) {
                    continue;
                }
                try {
#ifdef LEMON_BENCHMARK
                    auto start = std::chrono::high_resolution_clock::now();
#endif
                    auto traj = chemfiles::Trajectory::memory_reader(pair.second.data(), pair.second.size(), "MMTF/GZ");
                    auto entry = traj.read();
                    results[th].emplace_back(
                        worker(std::move(entry), pair.first));
#ifdef LEMON_BENCHMARK
                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration =
                        std::chrono::duration_cast<std::chrono::microseconds>(
                            stop - start);
                    std::cerr << it->string() + "\t" + pair.first + "\t" +
                                     std::to_string(duration.count()) + "\n";
#endif
                } catch (...) {
                }
            }
        }
    };

    for (auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread(call_function, work_iter, work_iter + grainsize);
        work_iter += grainsize;
    }
    threads.back() = std::thread(call_function, work_iter, pathvec.end());

    for (auto&& i : threads) {
        i.join();
    }
    for (const auto& thread_result : results) {
        for (const auto& sub_result : thread_result.second) {
            collector(sub_result);
        }
    }
}

#else

template <typename Function, typename Collector>
inline void run_parallel(Function&& worker, const std::string& p,
                         Collector& collector, size_t ncpu = 1,
                         const Entries& entries = Entries(),
                         const Entries& skip_entries = Entries()) {
    using ret = typename std::result_of<Function&(chemfiles::Frame,
                                                  const std::string&)>::type;
    auto pathvec = read_hadoop_dir(p);
    thread_pool threads(ncpu);
    threaded_queue<std::list<ret>> results;

    for (const auto& path : pathvec) {
        threads.queue_task([path, &results, &worker, &entries, &skip_entries] {
            std::ifstream data(path, std::istream::binary);
            Hadoop sequence(data);
            std::list<ret> mini_collector;

            while (sequence.has_next()) {
                auto pair = sequence.next();
                if (entries.size() && entries.count(pair.first) == 0) {
                    continue;
                }
                if (skip_entries.size() &&
                    skip_entries.count(pair.first) != 0) {
                    continue;
                }
                try {
#ifdef LEMON_BENCHMARK
                    auto start = std::chrono::high_resolution_clock::now();
#endif
                    auto traj = chemfiles::Trajectory::memory_reader(pair.second.data(), pair.second.size(), "MMTF/GZ");
                    auto entry = traj.read();
                    mini_collector.emplace_back(
                        worker(std::move(entry), pair.first));
#ifdef LEMON_BENCHMARK
                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration =
                        std::chrono::duration_cast<std::chrono::microseconds>(
                            stop - start);
                    std::cerr << path.string() + "\t" + pair.first + "\t" +
                                     std::to_string(duration.count()) + "\n";
#endif
                } catch (...) {
                }
            }

            results.push_back(std::move(mini_collector));
        });
    }

    std::size_t tasks_complete = 0;
    while (auto result = results.pop_front()) {
        for (auto sub_result : *result)
            collector(sub_result);
        ++tasks_complete;
        if (tasks_complete == pathvec.size()) {
            break;
        }
    }
}

#endif // LEMON_USE_ASYNC

} // namespace lemon
#endif
