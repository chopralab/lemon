#include "chemfiles.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-prototypes"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-braces"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wsign-conversion"

#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;

template<typename T>
bool check(const chemfiles::optional<T>& o) {
    return o != chemfiles::nullopt;
}

template<typename T>
T get(const chemfiles::optional<T>& o) {
    if (o) {
        return *o;
    } else {
        throw pybind11::index_error("Cannot dereference nullopt");
    }
    chemfiles::unreachable();
}

template<typename T, typename ret, int maxval>
ret get_index(const T& v, int i) {
    if (i >= 0) {
        if (i >= maxval) {
            throw pybind11::index_error("index is too large");
        }
        return v[static_cast<size_t>(i)];
    }
    if (-i > maxval) {
        throw pybind11::index_error("index is too small");
    }
    return v[static_cast<size_t>(maxval + i)];
}

template<typename T, typename ret>
ret get_index(const T& v, int i) {
    if (i >= 0) {
        if (i >= static_cast<int>(v.size())) {
            throw pybind11::index_error("index is too large");
        }
        return v[static_cast<size_t>(i)];
    }
    if (-i > static_cast<int>(v.size())) {
        throw pybind11::index_error("index is too small");
    }
    return v[static_cast<size_t>(static_cast<int>(v.size()) + i)];
}

chemfiles::Frame open_model_in_file(const std::string& filename, size_t index) {
    chemfiles::Trajectory traj(filename);
    return traj.read_step(index);
}

chemfiles::Frame* open_file(const std::string& filename) {
    chemfiles::Trajectory traj(filename);
    return new chemfiles::Frame(traj.read());
}

void write_file(const chemfiles::Frame& frame, const std::string& filename) {
    chemfiles::Trajectory traj(filename, 'w');
    traj.write(frame);
}

void append_file(const chemfiles::Frame& frame, const std::string& filename) {
    chemfiles::Trajectory traj(filename, 'a');
    traj.write(frame);
}

void add_chemfiles_features(py::module& m) {
    using namespace chemfiles;

    /**************************************************************************
     * Vector-likes
     **************************************************************************/
    py::class_<Bond>(m, "Bond")
        .def("__getitem__", &get_index<Bond,size_t,2>);

    py::bind_vector<std::vector<Bond>>(m, "BondVec");

    py::class_<Angle>(m, "Angle")
        .def("__getitem__", &get_index<Angle,size_t,3>);

    py::bind_vector<std::vector<Angle>>(m, "AngleVec");

    py::class_<Dihedral>(m, "Dihedral")
        .def("__getitem__", &get_index<Dihedral,size_t,4>);

    py::bind_vector<std::vector<Dihedral>>(m, "DihedralVec");

    py::class_<Improper>(m, "Improper")
        .def("__getitem__", &get_index<Improper,size_t,4>);

    py::bind_vector<std::vector<Improper>>(m, "ImproperVec");

    py::class_<Vector3D>(m, "Vector3D")
        .def(py::init<double,double,double>())
        .def("__getitem__", &get_index<Vector3D,double,3>)
        .def("norm", &Vector3D::norm);

    py::class_<Matrix3D>(m, "Matrix3D")
        .def(py::init<double,double,double,double,double,double,double,double,double>())
        .def("__getitem__", [](const Matrix3D& m, int i){
            if (static_cast<size_t>(i) >= 3) {
                throw pybind11::index_error("invalid matrix index");
            }

            return m[static_cast<size_t>(i)];
        });

    py::bind_vector<std::vector<Vector3D>>(m, "Coordinates");

    py::class_<span<Vector3D>>(m, "CoordinateSpan");

    /**************************************************************************
     * Optional
     **************************************************************************/
    py::class_<optional<uint64_t>> (m, "OptionalUInt64")
        .def("check", check<uint64_t>)
        .def("get", get<uint64_t>);

    py::class_<optional<std::string>> (m, "OptionalString")
        .def("check", check<std::string>)
        .def("get", get<std::string>);

    py::class_<optional<double>> (m, "OptionalDouble")
        .def("check", check<double>)
        .def("get", get<double>);

    py::class_<optional<const Property&>> (m, "OptionalProperty")
        .def("check", check<const Property&>)
        .def("get", get<const Property&>,
            py::return_value_policy::reference_internal);

    py::class_<optional<const Residue&>> (m, "OptionalResidue")
        .def("check", check<const Residue&>)
        .def("get", get<const Residue&>,
            py::return_value_policy::reference_internal);

    /**************************************************************************
     * Property
     **************************************************************************/
    py::class_<Property> (m, "Property")
        .def(py::init<bool>())
        .def(py::init<std::string>())
        .def(py::init<Vector3D>())
        .def(py::init<int>())
        .def(py::init<double>())
        .def("kind", &Property::kind)
        .def("as_bool", &Property::as_bool)
        .def("as_double", &Property::as_double)
        .def("as_vector3d", &Property::as_vector3d)
        .def("as_string", &Property::as_string);

    py::enum_<Property::Kind>(m, "Kind")
        .value("BOOL", Property::BOOL)
        .value("DOUBLE", Property::DOUBLE)
        .value("STRING", Property::STRING)
        .value("VECTOR3D", Property::VECTOR3D);

    /**************************************************************************
     * Residue
     **************************************************************************/
    chemfiles::optional<const chemfiles::Property&> (Residue::*residue_get)
        (const std::string&) const = &Residue::get;

    py::class_<Residue>(m, "Residue")
        .def(py::init<std::string>())
        .def(py::init<std::string, size_t>())
        .def("__len__", &Residue::size)
        .def("__iter__", [](const Residue& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("name", &Residue::name)
        .def("contains", &Residue::contains)
        .def("id", &Residue::id)
        .def("get", residue_get);

    py::bind_vector<std::vector<Residue>>(m, "ResidueVec");

    /**************************************************************************
     * Topology
     **************************************************************************/
    py::class_<Topology>(m, "Topology")
        .def("__len__", &Topology::size)
        .def("__get_item__", get_index<Topology,const Atom&>,
            py::return_value_policy::reference_internal)
        .def("__iter__", [](const Topology& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("residue", &Topology::residue,
            py::return_value_policy::reference_internal)
        .def("residues", &Topology::residues,
            py::return_value_policy::reference_internal)
        .def("residue_for_atom", &Topology::residue_for_atom)
        .def("are_linked", &Topology::are_linked)
        .def("bonds", &Topology::bonds,
            py::return_value_policy::reference_internal)
        .def("angles", &Topology::angles,
            py::return_value_policy::reference_internal)
        .def("dihedrals", &Topology::dihedrals,
            py::return_value_policy::reference_internal)
        .def("impropers", &Topology::impropers,
            py::return_value_policy::reference_internal);

    /**************************************************************************
     * Frame
     **************************************************************************/
    chemfiles::optional<const chemfiles::Property&> (Frame::*frame_get)
        (const std::string&) const = &Frame::get;

    py::class_<Frame>(m, "Frame")
        .def(py::init<>())
        .def("__len__", &Frame::size)
        .def("__get_item__", get_index<Frame,const Atom&>,
            py::return_value_policy::reference_internal)
        .def("__iter__", [](const Frame& v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("positions", [](Frame& v) {
            return v.positions();
        }, py::keep_alive<0, 1>())
        .def("positions_const", [](const Frame& v) {
            return v.positions();
        },  py::return_value_policy::reference_internal)
        .def("topology", &Frame::topology,
            py::return_value_policy::reference_internal)
        .def("distance", &Frame::distance)
        .def("angle", &Frame::angle)
        .def("dihedral", &Frame::dihedral)
        .def("out_of_plane", &Frame::out_of_plane)
        .def("get", frame_get);

    /**************************************************************************
     * Atom
     **************************************************************************/
    chemfiles::optional<const chemfiles::Property&> (Atom::*atom_get)
        (const std::string&) const = &Atom::get;

    py::class_<Atom>(m, "Atom")
        .def("name", &Atom::name)
        .def("type", &Atom::type)
        .def("mass", &Atom::mass)
        .def("charge", &Atom::charge)
        .def("full_name", &Atom::full_name)
        .def("vdw_radius", &Atom::vdw_radius)
        .def("covalent_radius", &Atom::covalent_radius)
        .def("atomic_number", &Atom::atomic_number)
        .def("get", atom_get);

    /**************************************************************************
    * File IO
    ***************************************************************************/
    m.def("open_model_in_file", open_model_in_file);
    m.def("open_file", open_file);
    m.def("write_file", write_file);
    m.def("append_file", append_file);
}

#pragma clang diagnostic pop
