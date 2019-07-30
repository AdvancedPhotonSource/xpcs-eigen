#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/numpy.hpp>

#include <xpcs/configuration.h>

namespace np = boost::python::numpy;
namespace python = boost::python;

template<typename T>
T * getNPData(python::object np_array){
    np::ndarray a = extract<np::ndarray>(np_array);
    return T * a.get_data();
}

void Configuration::from_pydict(python::dict &xpcs) {
    dqmap = getNPData<int>(xpcs.get("dqmap"));
    sqmap = getNPData<int>(xpcs.get("sqmap"));

    frameStart = extract<int>(xpcs.get("data_begin"));
    frameEnd = extract<int>(xpcs.get("data_end"));
    frameStartTodo = extract<int>(xpcs.get("data_begin_todo"));
    frameEndTodo = extract<int>(xpcs.get("data_end_todo"));
    delays_per_level_ = extract<int>(xpcs.get("delays_per_level"));
    darkFrameStart = extract<int>(xpcs.get("dark_begin_todo"));
    darkFrameEnd = extract<int>(xpcs.get("dark_end_todo"));

    two2one_window_size_ =
        extract<int>(xpcs.get("twotime2onetime_window_size"));
    if (two2one_window_size_ == NULL) two2one_window_size_ = 1;

    frame_stride_ = extract<long>(xpcs.get("stride_frames"));
    frame_average_ = extract<long>(xpcs.get("avg_frames"));

    normalizedByFrameSum = extract<bool>(xpcs.get("normalize_by_framesum"));
    if (normalizedByFrameSum == NULL) normalizedByFrameSum = false;

    if (darkFrameStart == darkFrameEnd || darkFrameEnd == 0) {
        darkFrameStart = 0;
        darkFrameEnd = 0;
        darkFrames = 0;
    } else {
        darkFrameStart = darkFrameStart - 1;
        darkFrameEnd = darkFrameEnd - 1;
        darkFrames = darkFrameEnd - darkFrameStart + 1;

        darkThreshold = extract<float>(xpcs.get("lld"));
        darkSigma = extract<float>(xpcs.get("sigma"));
    }

    m_totalStaticPartitions = 0;
    m_totalDynamicPartitions = 0;

    // detector shit
    python::dict detector = extract<python::dict>(xpcs.get("detector"));
    xdim = extract<int>(detector.get("x_dimension"));
    ydim = extract<int>(detector.get("y_dimension"));
    m_detDpixX = extract<float>(detector.get("x_pixel_size"));
    m_detDpixY = extract<float>(detector.get("y_pixel_size"));
    m_detAdhupPhot = extract<float>(detector.get("adu_per_photon"));
    m_detPreset = extract<float>(detector.get("exposure_time"));
    m_detEfficiency = extract<float>(detector.get("efficiency"));

    float det_dist = extract<float>(detector.get("distance"));
    float flux = extract<float>(detector.get("beam_intensity_transmitted"));
    float thickness = extract<float>(xpcs.get("thickness"));

    m_normFactor = 1. / m_detEfficiency / m_detAdhupPhot / m_detPreset;
    m_normFactor /= (m_detDpixX / det_dist * m_detDpixY / det_dist);

    m_staticWindow = extract<int>(xpcs.get("static_mean_window_size"));
    flatfieldEnabled = extract<bool>(xpcs.get("flatfield_enabled"));
    if (flatfield_enabled != NULL) {
        flatfield = static_cast<double *>(
            extract<np::ndarray>(detector.get("flatfield")).get_data());
    } else {
        flatfield = new double[xdim * ydim];
        for (int i = 0; i < xdim * ydim; i++) flatfield[i] = 1.;
    }

    twotime_ = false; // TODO move default values to constructor
    std::string str = extract<std::string>(xpcs.get("analysis_type"));
    if (str.compare("Twotime") == 0) {
        twotime_ = true;
        np::ndarray qphi_bins = extract<np::ndarray>(xpcs.get("qphi_bin_to_process"));
        int size = qphi_bins.shape(0);
        int * qphibins = reinterpret_cast<int *>(qphi_bins.get_data());
        for (int i = 0; i < size; i++) 
            qphi_bin_to_process_.insert(qphi_bin_to_process_.begin(), qphibins, qphibins + size);
    } 
}
