#include <boost/python.hpp>

using namespace boost::pyhton;

#include <xpcs/configuration.h>

BOOST_PYTHON_MODULE(xpcs_eigen) {
    class<Configuration>("Configuration")
        .def("init", &Configuration::from_pydict())
        ;
}
