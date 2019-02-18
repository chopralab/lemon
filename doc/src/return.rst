.. _return:

Handling Custom Returns in C++
==============================

The **C++** API allows users to return custom objects from their *workflows*.

To do so, the user must create a *functor* object which overloads the
`operator()` function. Lemon provides two implementations for **C++** streams
(such as `std::cout`) and associative arrays (such as `std::map`). These two
examples are given below to guide users in the creation of their own custom
functors for new types. It is recommended that the user use **C++** templates
to obtain the most flexible implementations.

.. doxygenstruct:: lemon::print_combine
    :members:

.. doxygenstruct:: lemon::map_combine
    :members:

Python alternative
==================

This feature is under current development for the **Python** API, but is not
yet supported. Instead, users should store the results of their workflows in
the `self` object and rely on the `finalize` function to analyze the result.
