#include "lemon/options.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Default Options") {
    const char* argv[] = {"junk"};

	lemon::Options opts(1, argv);
    CHECK(opts.work_dir() == ".");
    CHECK(opts.ncpu() == 1);
    CHECK(opts.entries().empty());
    CHECK(opts.skip_entries().empty());
}
