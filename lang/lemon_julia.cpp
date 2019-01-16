#include <julia.h>
#include <string>
#include <iostream>

#include "lemon/lemon.hpp"

//JULIA_DEFINE_FAST_TLS()

int main(int argc, char *argv[]) {
    lemon::Options o;
    std::string julia_script("lemon.jl");
    o.add_option("julia_script,j", julia_script, "Julia script to run");
    o.parse_command_line(argc, argv);

    /* required: setup the Julia context */
    jl_init();

    /* run Julia commands */
    std::string julia_command("include(\"");
    julia_command += julia_script;
    julia_command += "\")";
    jl_eval_string("@eval Main import Base.MainInclude: eval, include");
    jl_eval_string(julia_command.c_str());

    if (jl_exception_occurred()) {
        std::cerr << "Error in provided script: "
                  << jl_typeof_str(jl_exception_occurred()) << std::endl;
        return 1;
    }

    auto mod = jl_eval_string("Lemon");

    if (!mod || jl_exception_occurred()) {
        std::cerr << "Error in getting Lemon module: "
                  << jl_typeof_str(jl_exception_occurred()) << std::endl;
        return 1;
    }

    auto func = jl_get_function(reinterpret_cast<jl_module_t*>(mod),"worker");

    if (!func) {
        std::cerr << "Error loading worker function." << std::endl;
        return 2;
    }

    auto worker = [&func](chemfiles::Frame complex, const std::string&) {
        auto res0 = jl_call0(func);
        return std::string("");
        ///auto jl_frame = jl_box_voidpointer(asdf2);
        //auto res = jl_call1(func, jl_frame);

        //return std::string("");
    };

    //lemon::launch<lemon::print_combine>(o, worker, std::cout);

    jl_atexit_hook(0);
    return 0;
}
