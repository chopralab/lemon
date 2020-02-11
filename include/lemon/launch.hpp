#ifndef LEMON_LAUNCH_HPP
#define LEMON_LAUNCH_HPP

#include "lemon/constants.hpp"
#include "lemon/options.hpp"
#include "lemon/parallel.hpp"

#include <iostream>
#include <ostream>

namespace lemon {

class LemonPythonBase {
  public:
    virtual ~LemonPythonBase() = default;
    
    LemonPythonBase() = default;
    LemonPythonBase(const LemonPythonBase& other) = delete;
    LemonPythonBase(LemonPythonBase&& other) noexcept = delete;

    LemonPythonBase& operator=(const LemonPythonBase& other) = delete;
    LemonPythonBase& operator=(LemonPythonBase&& other) noexcept = delete;


    virtual std::string worker(const chemfiles::Frame*, const std::string&) = 0;
    virtual void finalize() { /*Do nothing*/
    }
};

//! Functor to combine the results of a workflow into a map
//!
//! This is a template functor which opens combines the results of a workflow
//! into a single map object. Use this class only when the work flow returns a
//! map oject or other type of associative array.
template <typename Map1> struct map_combine {

    //! Construct a map_combine class which fills the `collector`
    //!
    //! This constructor will copy all results to the `collector` object. This
    //! operation is performed when the workflow completes for all entries.
    //! \param collector A map like object to store results in.
    map_combine(Map1& collector) : internal_map_(collector) {}

    //! Add a map to the current `collector`.
    template <typename Map2 = Map1> void operator()(const Map2& map2) const {
        for (const auto& sc : map2) {
            internal_map_[sc.first] += sc.second;
        }
    }

  private:
    Map1& internal_map_;
};

//! Functor to stream workflow results to an ostream object
//!
//! This is a template functor which streams workflow results to a `collector`.
//! This class is typically used to stream text results from workflows to the
//! `std::cout` object. Any class which derives from `std::ostream` is
//! supported.
struct print_combine {

    //! Construct a print_combine which streams results to `collector`.
    //!
    //! This contructor will create a `print_combine` object which streams
    //! results to `collector`. This argument must derive from the
    //! `std::ostream` class and therefore provide an overload to the `<<`
    //! operator. Common examples of `collector` are `std::cout` and friends.
    //! All workflow results are streamed to the `collect` after all entries
    //! have been evaluated. \param collector An object for each workflow
    //! objects will be streamed to
    print_combine(std::ostream& collector) : internal_stream_(collector) {}

    template <typename Res> void operator()(const Res& map) const {
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
template <typename Function, typename Collector>
int launch(const Options& o, Function&& worker, Collector& collect) {
    const auto& p = o.work_dir();
    auto threads = o.ncpu();
    const auto& entries = read_entry_file(o.entries());
    const auto& skip_entries = read_entry_file(o.skip_entries());

    try {
        lemon::run_parallel(worker, p, collect, threads, entries, skip_entries);
    } catch (std::runtime_error& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}

} // namespace lemon

#endif
