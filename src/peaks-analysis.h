/*
    Copyright (C) CENATAV, DATYS - All Rights Reserved
    Unauthorized copying of this file, via any medium is strictly prohibited
    Proprietary and confidential
    Written by Victor M. Mendiola Lau <vmendiola@cenatav.co.cu>, March 2017
*/

#ifndef PEAKS_VALLEYS_H_INCLUDED
#define PEAKS_VALLEYS_H_INCLUDED

#include "config.hpp"

/// <summary>
/// Represents the monotony of the function in interval I
/// </summary>
enum Monotony
{
    GROW,       // the function increases
    ABATE,      // the function decreases
    STABLE      // the function is stable
};

/// <summary>
/// Represents the possible extreme points
/// </summary>
enum ExtremeType
{
    MIN,        // the extreme point is a MIN
    MAX,        // the extreme point is a MAX
    NONE        // no extreme point
};

/// <summary>
/// Holds peak related information
/// </summary>
struct PeakInfo
{
    int lower_bound;

    int upper_bound;

    int height_index;

    double peak_height;

    double peak_area;

    PeakInfo()
    {
        PeakInfo(0, 0, 0, 0.0, 0.0);
    }

    PeakInfo(int lb, int ub, int hi, double ph, double pa)
    {
        // initializing struct fields
        lower_bound = lb;
        upper_bound = ub;
        height_index = hi;
        peak_height = ph;
        peak_area = pa;
    }
};

// ----------------------------------------------------

/// <summary>
/// Computes the slope of for the given angle
/// </summary>
double compute_slope(double angle);

// ----------------------------------------------------

/// <summary>
/// Finds peaks and valleys in certain histogram
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <returns>Histogram's peaks and valleys</returns>
/// <remarks>Does not perform histogram normalization</remarks>
std::vector<PeakInfo> find_in_histogram(std::vector<double> histogram, int dx, int smoothness, double growth_angle, double abate_angle, double height_thres);

/// <summary>
/// Classifies an extreme point according to monotony
/// </summary>
/// <param name="previous">Monotony of the previous interval</param>
/// <param name="current">Monotony of the current interval</param>
/// <returns>Extreme point type</returns>
ExtremeType get_extreme_type(Monotony previous, Monotony current);

/// <summary>
/// Finds the monotony in an interval using central differences
/// </summary>
/// <param name="fa">Function value at a [f(a)]</param>
/// <param name="fb">Function value at b [f(b)]</param>
/// <param name="dx">Interval size</param>
/// <param name="grow_thres">Threshold for a growing function</param>
/// <param name="abate_thres">Threshold for an abating function</param>
/// <returns>The type of monotony</returns>
Monotony find_monotony(double fa, double fb, int dx, double growth_thres, double abate_thres);

/// <summary>
/// Computes peak height
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <param name="peak_max_idx">Index at which peak reaches its max value</param>
/// <returns>Peak height</returns>
/// <remarks>Does not perform histogram normalization</remarks>
double compute_peak_height(std::vector<double> histogram, int lb, int ub, int *peak_max_idx);

/// <summary>
/// Computes peak area
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <returns>Peak area computed</returns>
/// <remarks>Does not perform histogram normalization</remarks>
double compute_peak_area(std::vector<double> histogram, int lb, int ub);

// area of type A1
double compute_peak_area1(std::vector<double> histogram, int lb, int ub);

// area of type A2
double compute_peak_area2(std::vector<double> histogram, int lb, int ub);

/// <summary>
/// Computes peak statistics
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <returns>Peak statistics</returns>
PeakInfo compute_peak_statistics(std::vector<double> histogram, int lb, int ub);

// ---------------------------------- C-API members ----------------------------------

/// <summary>
/// Finds peaks and valleys in certain histogram
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="h_length">Histogram length</param>
/// <param name="dx">Interval size (dx) (size of interval of analysis I)</param>
/// <param name="smoothness">Indicates the "smoothness" of the curve</param>
/// <param name="growth_angle">Growth angle of certain interval (indicates if the function is really growing)</param>
/// <param name="abate_angle">Abate angle of certain interval (indicates if the function is really abating)</param>
/// <param name="height_thres">Threshold for filtering peaks that do not match required height</param>
/// <param name="peaks_count">Amount of peaks found</param>
/// <param name="valleys_count">Amount of valleys found</param>
/// <returns>Histogram's peaks and valleys</returns>
/// <remarks>Does not perform histogram normalization</remarks>
API_FUNC(PeakInfo*) find_in_histogram(double* histogram, int h_length, int dx, int smoothness, double growth_angle, double abate_angle, double height_thres, int* peaks_count);

/// <summary>
/// Computes peak statistics for further analysis
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <returns>Peak height, area, etc.</returns>
/// <remarks>Does not perform histogram normalization</remarks>
API_FUNC(PeakInfo) compute_peak_statistics(double* histogram, int lb, int ub);

/// <summary>
/// Frees a 'dangling' pointer of PeakInfo
/// </summary>
/// <param name="ptr">The pointer in question</param>
/// <remarks>This should be called after calling C-API 'find_in_histogram' function</remarks>
API_FUNC(void) delete_peak_info_ptr(PeakInfo* ptr);

#endif // PEAKS_VALLEYS_H_INCLUDED