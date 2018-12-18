#include "lemon/options.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Default Options") {
    const char* argv[] = {"junk"};

	lemon::Options opts(1, argv);
    CHECK(opts.work_dir() == ".");
    CHECK(opts.ncpu() == 1);
    CHECK(opts.entries().empty());

	std::string junk;
    CHECK_THROWS(opts.add_option("junk", junk));
    CHECK_THROWS(opts.parse_command_line(1, argv));
}

TEST_CASE("Custom Options") {
    const char* argv[] = {"junk","-w","home","-d","7.0","-n","2","-e","entries"};
    int argc = 9;

    double distance = 0;
    std::string reference;
    lemon::Options opts;
    opts.add_option("distance,d", distance);
    opts.add_option("reference,r", reference);

    opts.parse_command_line(argc, argv);
    CHECK(opts.work_dir() == "home");
    CHECK(distance == 7.0);
    CHECK(reference.empty());
    CHECK(opts.ncpu() == 2);
    CHECK(opts.entries() == "entries");

	CHECK_THROWS(opts.add_option("junk", reference));
    CHECK_THROWS(opts.parse_command_line(argc, argv));
}
