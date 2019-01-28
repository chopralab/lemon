#ifndef LEMON_LAUNCH_HPP
#define LEMON_LAUNCH_HPP

#include "lemon/options.hpp"
#include "lemon/parallel.hpp"
#include "lemon/constants.hpp"

#include <ostream>

namespace lemon {

class LemonPythonBase : public boost::noncopyable {
public:
    virtual ~LemonPythonBase() {}
    virtual std::string worker(const chemfiles::Frame*, const std::string&) = 0;
    virtual void finalize() {/*Do nothing*/}
};

template<typename Map1>
struct map_combine {

    map_combine(Map1& collector):
        internal_map_(collector) {}

    template<typename Map2 = Map1>
    void operator()(const Map2& map2) const {
        for (const auto& sc : map2) {
            internal_map_[sc.first] += sc.second;
        }
    }

private:
    Map1& internal_map_;
};

struct print_combine {

    print_combine(std::ostream& collector):
        internal_stream_(collector){}

    template<typename Res>
    void operator()(const Res& map) const {
        internal_stream_ << map;
    }

private:
    std::ostream& internal_stream_;
};

//! Launch a **Lemon** workflow.
//!
//! This function reads **Lemon** options and passes them to the appropriate
//! `run_parallel` function. This is the main entry point of a C++ **Lemon**
//! program.
//! \param [in] o An instance of the `Options` used to pass arguments to Lemon
//! \param worker Function object representing the body of the workflow.
//! \param collect Function object for collect the results of `worker`.
//! \return 0 on success or a non-zero integer on error.
template<typename Function, typename Collector>
int launch(const Options& o, Function&& worker, Collector& collect) {
    auto p = o.work_dir();
    auto threads = o.ncpu();
    auto entries = read_entry_file(o.entries());
    auto skip_entries = read_entry_file(o.skip_entries());

    try {
        lemon::run_parallel(worker, p, collect, threads, entries, skip_entries);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}

} // namespace lemon

#endif
