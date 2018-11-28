#ifndef LEMON_PARALLEL_HPP
#define LEMON_PARALLEL_HPP

#include "lemon/hadoop.hpp"

#ifdef LEMON_USE_ASYNC
#include "lemon/thread_pool.hpp"
#endif

namespace lemon {

#ifndef LEMON_USE_ASYNC

// Old threading behavior
template <class Function>
inline void run_parallel(Function&& worker, const fs::path& p, size_t ncpu) {
    auto pathvec = read_hadoop_dir(p);
    std::vector<std::thread> threads(ncpu);

    // Total number of jobs for each thread
    const size_t grainsize = pathvec.size() / ncpu;
    auto work_iter = pathvec.begin();
    using iter = std::vector<fs::path>::iterator;

    auto call_function = [&worker](iter first, iter last) {
        for (auto it = first; it != last; ++it) {
            std::ifstream data(it->string(), std::istream::binary);
            Hadoop sequence(data);

            while (sequence.has_next()) {
                auto pair = sequence.next();
                const auto entry = std::string(pair.first.data() + 1, 4);
                try {
                    auto traj = chemfiles::Trajectory(std::move(pair.second),
                                                      "MMTF/GZ");
                    auto complex = traj.read();
                    worker(std::move(complex), entry);
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
}

template <typename Function, typename Combiner, typename ret>
inline void run_parallel(Function&& worker, Combiner&& combine,
                         const fs::path& p, ret& collector, size_t ncpu) {
    auto pathvec = read_hadoop_dir(p);
    std::vector<std::thread> threads(ncpu);

    // Total number of jobs for each thread
    const size_t grainsize = pathvec.size() / ncpu;
    auto work_iter = pathvec.begin();
    using iter = std::vector<fs::path>::iterator;
    std::unordered_map<std::thread::id, ret> results;

    auto call_function = [&worker, &results, &combine](iter first, iter last) {
        auto th = std::this_thread::get_id();
        for (auto it = first; it != last; ++it) {
            std::ifstream data(it->string(), std::istream::binary);
            Hadoop sequence(data);

            while (sequence.has_next()) {
                auto pair = sequence.next();
                const auto entry = std::string(pair.first.data() + 1, 4);
                try {
                    auto traj = chemfiles::Trajectory(std::move(pair.second),
                                                      "MMTF/GZ");
                    auto complex = traj.read();
                    combine(results[th], worker(std::move(complex), entry));
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
        combine(collector, thread_result.second);
    }
}

#else

template <class Function>
inline void run_parallel(Function&& worker, const fs::path& p, size_t ncpu) {
    auto pathvec = read_hadoop_dir(p);
    thread_pool threads(ncpu);
    threaded_queue<size_t> results;

    for (const auto& path : pathvec) {
        threads.queue_task([path, &results, &worker] {
            std::ifstream data(path.string(), std::istream::binary);
            Hadoop sequence(data);

            while (sequence.has_next()) {
                auto pair = sequence.next();
                const auto entry = std::string(pair.first.data() + 1, 4);
                try {
                    auto traj = chemfiles::Trajectory(std::move(pair.second),
                                                      "MMTF/GZ");
                    auto complex = traj.read();
                    worker(std::move(complex), entry);
                } catch (...) {
                }
            }

            results.push_back(0);
        });
    }

    std::size_t tasks_complete = 0;
    while (auto result = results.pop_front()) {
        ++tasks_complete;
        if (tasks_complete == pathvec.size()) {
            break;
        }
    }
}

template <typename Function, typename Combiner,
          typename ret = std::result_of_t<Function&()> >
inline void run_parallel(Function&& worker, Combiner&& combine,
                         const fs::path& p, ret& collector, size_t ncpu) {
    auto pathvec = read_hadoop_dir(p);
    thread_pool threads(ncpu);
    threaded_queue<ret> results;

    for (const auto& path : pathvec) {
        threads.queue_task([path, &results, &worker, &combine] {
            std::ifstream data(path.string(), std::istream::binary);
            Hadoop sequence(data);
            ret mini_collector;

            while (sequence.has_next()) {
                auto pair = sequence.next();
                const auto entry = std::string(pair.first.data() + 1, 4);
                try {
                    auto traj = chemfiles::Trajectory(std::move(pair.second),
                                                      "MMTF/GZ");
                    auto complex = traj.read();
                    combine(mini_collector, worker(std::move(complex), entry));
                } catch (...) {
                }
            }

            results.push_back(mini_collector);
        });
    }

    std::size_t tasks_complete = 0;
    while (auto result = results.pop_front()) {
        combine(collector, *result);
        ++tasks_complete;
        if (tasks_complete == pathvec.size()) {
            break;
        }
    }
}

#endif  // LEMON_USE_ASYNC

}  // namespace lemon
#endif
