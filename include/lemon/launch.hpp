#ifndef LEMON_LAUNCH_HPP
#define LEMON_LAUNCH_HPP

#include "lemon/options.hpp"
#include "lemon/parallel.hpp"
#include "lemon/constants.hpp"

namespace lemon {

template<typename Map1, typename Map2>
struct map_combine {
    void operator()(Map1& map1, const Map2& map2) const {
        for (const auto& sc : map2) {
            map1[sc.first] += sc.second;
        }
    }
};

template<typename OS, typename Map>
struct print_combine {
    void operator()(OS& os, const Map& map) const {
        os << map;
    }
};

template< template<typename C, typename R> class Combiner, typename Function,
         typename Collector>
int launch(const Options& o, Function&& worker, Collector& collect) {
    auto p = o.work_dir();
    auto threads = o.ncpu();
    auto entries = read_entry_file(o.entries());
    auto se = o.skip_entries();

    std::unordered_set<std::string> skip_entries;
    if (se.empty()) {
        skip_entries = large_entries;
    } else if (se != "*none*") {
        skip_entries = read_entry_file(se);
    }

    using ret = typename std::result_of<Function&(chemfiles::Frame,
                                                  const std::string&)>::type;
    Combiner<Collector, ret> combiner;

    try {
        lemon::run_parallel(worker, combiner, p, collect, threads, entries, skip_entries);
    } catch(std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}

} // namespace lemon

#endif
