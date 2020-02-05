#ifndef LEMON_THREAD_POOL_HPP
#define LEMON_THREAD_POOL_HPP

#include <algorithm>
#include <chrono>
#include <deque>
#include <future>
#include <mutex>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>

#include "lemon/external/gaurd.hpp"

LEMON_EXTERNAL_FILE_PUSH
#include <chemfiles/external/optional.hpp>
LEMON_EXTERNAL_FILE_POP

namespace lemon {

template <class T> struct threaded_queue {
    using lock = std::unique_lock<std::mutex>;
    void push_back(T t) {
        {
            lock l(m);
            data.push_back(std::move(t));
        }
        cv.notify_one();
    }

    chemfiles::optional<T> pop_front() {
        lock l(m);
        cv.wait(l, [this] { return abort || !data.empty(); });
        if (abort)
            return chemfiles::nullopt;
        auto r = std::move(data.back());
        data.pop_back();
        return std::move(r);
    }

    void terminate() {
        {
            lock l(m);
            abort = true;
            data.clear();
        }
        cv.notify_all();
    }

    ~threaded_queue() { terminate(); }

  private:
    std::mutex m;
    std::deque<T> data;
    std::condition_variable cv;
    bool abort = false;
};

struct thread_pool {
    thread_pool(std::size_t n = 1) { start_thread(n); }
    thread_pool(thread_pool&&) = delete;
    thread_pool& operator=(thread_pool&&) = delete;
    ~thread_pool() = default;

    template <class F, class R = std::result_of_t<F&()>>
    std::future<R> queue_task(F task) {
        std::packaged_task<R()> p(std::move(task));
        auto r = p.get_future();
        tasks.push_back(std::move(p));
        return r;
    }

    template <class F, class R = std::result_of_t<F&()>>
    std::future<R> run_task(F task) {
        if (threads_active() >= total_threads()) {
            start_thread();
        }
        return queue_task(std::move(task));
    }

    void terminate() { tasks.terminate(); }

    std::size_t threads_active() const { return active; }

    std::size_t total_threads() const { return threads.size(); }

    void clear_threads() {
        terminate();
        threads.clear();
    }

    void start_thread(std::size_t n = 1) {
        while (n-- > 0) {
            threads.push_back(std::async(std::launch::async, [this] {
                while (auto task = tasks.pop_front()) {
                    ++active;
                    try {
                        (*task)();
                    } catch (...) {
                        --active;
                        throw;
                    }
                    --active;
                }
            }));
        }
    }

  private:
    std::vector<std::future<void>> threads;
    threaded_queue<std::packaged_task<void()>> tasks;
    std::atomic<std::size_t> active;
};
} // namespace lemon

#endif
