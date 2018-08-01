#include "lemon/options.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Options Test") {
    const char* argv[] = {"junk","-w","home","-d","7.0","-n","2","-e","entries"};
    int argc = 9;

    lemon::Options opts(argc, argv);
    CHECK(opts.work_dir() == "home");
    CHECK(opts.distance() == 7.0);
    CHECK(opts.reference().empty());
    CHECK(opts.npu() == 2);
    CHECK(opts.output() == ".");
    CHECK(opts.entries() == "entries");
}
