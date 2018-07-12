#ifndef RUN_HPP
#define RUN_HPP

#include <iterator>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

namespace benchmarker {
template <class Function, class iter>
void call_function(Function&& f, iter begin, iter end, size_t threadid) {
    for (auto it = begin; it != end; ++it) {
        const std::string pdbid = std::string(it->data(), 4);

        boost::filesystem::path entry(std::string(it->data(), 1));
        entry /= std::string(it->data() + 1, 1);
        entry /= pdbid + std::string(".mmtf.gz");

        try {
            auto traj = chemfiles::Trajectory(entry.string());

            if (traj.nsteps() > 1) {
                std::cerr << "File: " << entry
                          << " contrains more than more model. Only "
                             "considering the first"
                          << std::endl;
            }

            auto complex = traj.read();
            f(complex, pdbid, threadid);
        } catch (const std::range_error& e) {
            std::cerr << "Odd residue name in " << entry << std::endl;
        } catch (const std::length_error& e) {
            std::cerr << "Long residue name in " << entry << std::endl;
        } catch (const chemfiles::FormatError& e) {
            std::cerr << "Unsupported file format for " << entry << ". Skipping"
                      << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Unknown error: " << e.what() << " for " << entry
                      << std::endl;
        } catch (...) {
            std::cerr << "Unknown error for " << entry << std::endl;
        }
    }
}

template <class Function, class container>
void call_multithreaded(Function&& worker, const container& vec, size_t ncpu) {
    std::vector<std::thread> threads(ncpu);
    const int grainsize = vec.size() / ncpu;

    auto work_iter = std::cbegin(vec);
    size_t id = 0;
    for (auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread([&] {
            call_function(worker, work_iter, work_iter + grainsize, id);
        });
        ++id;
        work_iter += grainsize;
    }
    threads.back() = std::thread(
        [&] { call_function(worker, work_iter, work_iter + grainsize, id); });

    for (auto&& i : threads) {
        i.join();
    }
}
}

#endif
