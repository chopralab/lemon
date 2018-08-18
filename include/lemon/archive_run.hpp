#ifndef ARCHIVE_RUN_HPP
#define ARCHIVE_RUN_HPP

#include <iterator>
#include <string>
#include <vector>
#include <thread>

#include <boost/filesystem.hpp>

namespace lemon {
template <class Function, class iter>
void call_function(Function&& f, iter begin, iter end) {
    for (auto it = begin; it != end; ++it) {
        const std::string pdbid = std::string(it->data(), 4);

        boost::filesystem::path entry(std::string(it->data(), 1));
        entry /= std::string(it->data() + 1, 1);
        entry /= pdbid + std::string(".mmtf.gz");

        try {
            auto traj = chemfiles::Trajectory(entry.string());

#ifndef NDEBUG
            if (traj.nsteps() > 1) {
                std::cerr << "File: " << entry
                          << " contrains more than more model. Only "
                             "considering the first"
                          << std::endl;
            }
#endif

            auto complex = traj.read();
            f(std::move(complex), pdbid);
        } catch (...) {
            std::cerr << "Unknown error for " << entry << std::endl;
        }
    }
}

template <class Function, class container>
void run_archive(Function&& worker, const container& vec, size_t ncpu, size_t chunk = 1) {
    std::vector<std::thread> threads(ncpu);

    // Total number of jobs for each thread
    const size_t grainsize = vec.size() / ncpu;

    // Total number of jobs submitted to be run each time a thread is spun
    const size_t chunksize = grainsize / chunk;

    auto work_iter = vec.cbegin();

    for (size_t chunk_id = 1; chunk_id < chunk; ++chunk_id) {
        for (auto it = std::begin(threads); it != std::end(threads); ++it) {
            *it = std::thread([=] {
                call_function(worker, work_iter, work_iter + chunksize);
            });
            work_iter += chunksize;
        }

        for (auto&& i : threads) {
            i.join();
        }
    }

    for (auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread([=] {
            call_function(worker, work_iter, work_iter + chunksize);
        });
        work_iter += chunksize;
    }
    threads.back() = std::thread(
        [&] { call_function(worker, work_iter, vec.cend()); });

    for (auto&& i : threads) {
        i.join();
    }
}
}

#endif
