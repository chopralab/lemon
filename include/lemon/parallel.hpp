#ifndef LEMON_PARALLEL_HPP
#define LEMON_PARALLEL_HPP

#include "lemon/hadoop.hpp"

#ifdef LEMON_USE_ASYNC
#include "lemon/thread_pool.hpp"
#else
#include <thread>
#endif

namespace lemon {

#ifndef LEMON_USE_ASYNC

//! The `run_parallel` function launches jobs which do return data.
//!
//! Use this function to run the `worker` function on `ncpu` threads. The `worker`
//! should accept two arguments, a `chemfiles::Frame` and a `std::string`. It
//! must return a value as this value will be appended, using the `combine`
//! function object, the the `collector`.
//! See the `Lemon Workflow` documention for more details.
//! \param worker A function object (C++11 lambda, struct the with operator()
//!  overloaded, or std::function object) that the user wishes to apply.
//! \param combine A function object that combines the return value of `worker`
//!  to the `collector` object.
//! \param [in] p A path to the Hadoop sequence file directory.
//! \param [out] collector An object to hold the collection of values returned
//!  by the worker object which are appended with `combine`.
//! \param [in] ncpu The number of threads to use.
//! \param [in] entries Which entries to use. Not used if blank.
//! \param [in] skip_entries Which entries to skip. Not used if blank.
template <typename Function, typename Combiner, typename Collector>
inline void run_parallel(Function&& worker, Combiner&& combine,
                         const fs::path& p, Collector& collector, size_t ncpu = 1,
                         const Entries& entries = Entries(),
                         const Entries& skip_entries = Entries()) {
    auto pathvec = read_hadoop_dir(p);
    std::vector<std::thread> threads(ncpu);
    using ret = typename std::result_of<Function&(chemfiles::Frame, const std::string&)>::type;

    // Total number of jobs for each thread
    const int grainsize = static_cast<int>(pathvec.size() / ncpu);
    auto work_iter = pathvec.begin();
    using iter = std::vector<fs::path>::iterator;
    std::unordered_map<std::thread::id, std::list<ret>> results;

    auto call_function = [&worker, &results, &combine, &entries, &skip_entries](iter first, iter last) {
        auto th = std::this_thread::get_id();
        for (auto it = first; it != last; ++it) {
            std::ifstream data(it->string(), std::istream::binary);
            Hadoop sequence(data);

            while (sequence.has_next()) {
                auto pair = sequence.next();
                if (entries.size() && entries.count(pair.first) == 0) {
                    continue;
                }
                if (skip_entries.size() && skip_entries.count(pair.first) != 0) {
                    continue;
                }
                try {
                    auto traj = chemfiles::Trajectory(std::move(pair.second),
                                                      "MMTF/GZ");
                    auto complex = traj.read();
                    results[th].emplace_back(worker(std::move(complex), pair.first));
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
        for (auto sub_result : thread_result.second) {
            combine(collector, sub_result);
        }
    }
}

#else

template <typename Function, typename Combiner, typename Collector>
inline void run_parallel(Function&& worker, Combiner&& combine,
                         const fs::path& p, Collector& collector, size_t ncpu = 1,
                         const Entries& entries = Entries(),
                         const Entries& skip_entries = Entries()) {
    using ret = typename std::result_of<Function&(chemfiles::Frame, const std::string&)>::type;
    auto pathvec = read_hadoop_dir(p);
    thread_pool threads(ncpu);
    threaded_queue<std::list<ret>> results;

    for (const auto& path : pathvec) {
        threads.queue_task([path, &results, &worker, &combine, &entries, &skip_entries] {
            std::ifstream data(path.string(), std::istream::binary);
            Hadoop sequence(data);
            std::list<ret> mini_collector;

            while (sequence.has_next()) {
                auto pair = sequence.next();
                if (entries.size() && entries.count(pair.first) == 0) {
                    continue;
                }
                if (skip_entries.size() && skip_entries.count(pair.first) != 0) {
                    continue;
                }
                try {
                    auto traj = chemfiles::Trajectory(std::move(pair.second),
                                                      "MMTF/GZ");
                    auto complex = traj.read();
                    mini_collector.emplace_back(worker(std::move(complex), pair.first));
                } catch (...) {
                }
            }

            results.push_back(std::move(mini_collector));
        });
    }

    std::size_t tasks_complete = 0;
    while (auto result = results.pop_front()) {
        for (auto sub_result : *result)
            combine(collector, sub_result);
        ++tasks_complete;
        if (tasks_complete == pathvec.size()) {
            break;
        }
    }
}

#endif  // LEMON_USE_ASYNC

}  // namespace lemon
#endif
