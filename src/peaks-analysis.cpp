/*
    Copyright (C) CENATAV, DATYS - All Rights Reserved
    Unauthorized copying of this file, via any medium is strictly prohibited
    Proprietary and confidential
    Written by Victor M. Mendiola Lau <vmendiola@cenatav.co.cu>, March 2017
*/

#include <cmath>
#include <vector>
#include "peaks-analysis.h"

/// <summary>
/// Represents Math.PI
/// </summary>
const double pi = 3.141592653589793;

// ----------------------------------------------------

/// <summary>
/// Computes the slope of for the given angle
/// </summary>
double compute_slope(double angle)
{
    return tan(angle * pi / 180.0);
}

// ----------------------------------------------------

/// <summary>
/// Finds peaks and valleys in certain histogram
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <returns>Histogram's peaks and valleys</returns>
/// <remarks>Does not perform histogram normalization</remarks>
std::vector<PeakInfo> find_in_histogram(std::vector<double> histogram, int dx, int smoothness, double growth_angle, double abate_angle, double height_thres)
{
    // declaring variables that will hold the peaks
    std::vector<PeakInfo> peaks;

    // getting histogram length
    int histogram_length = (int)histogram.size();

    // this is because there's always a point shared between intervals
    int offset = dx - 1;

    // there can not exist peaks or valleys (there's at most one interval)
    if (dx < 2 || histogram_length < 2 * offset + 1)
        return peaks;

    // Ip: previous interval (Ip = [a, b])
    // Ic: current interval  (Ic = [b, c])
    int b = offset;
    int c = b + offset;

    // computing 'growth' and 'abate' slope thresholds
    double _growth_thres = compute_slope(growth_angle);
    double _abate_thres = compute_slope(abate_angle);

    // analyzing Ip monotony
    Monotony ip_monotony = find_monotony(histogram[0], histogram[b], dx, _growth_thres, _abate_thres);

    // setting previous extreme point
    ExtremeType previous_extreme = NONE;
    int previous_extreme_index = 1;

    if (ip_monotony == ABATE)
        previous_extreme = MAX;
    else if (ip_monotony == GROW)
        previous_extreme = MIN;

    // setting current extreme point
    ExtremeType current_extreme;
    int current_extreme_index = 0;

    // setting the first extreme point of the last registered shift
    ExtremeType previous_shift_extreme = NONE;
    int previous_shift_extreme_index = 0;

    // while the current interval is in range
    while (c < histogram_length)
    {
        // analyzing the monotony in Ic
        Monotony ic_monotony = find_monotony(histogram[b], histogram[c], dx, _growth_thres, _abate_thres);

        // if there was a change in the monotony => there's a extreme
        if (ip_monotony != ic_monotony)
        {
            // classifying the new extreme founded
            current_extreme = get_extreme_type(ip_monotony, ic_monotony);
            current_extreme_index = b;

            // if there was a shift indeed
            if (current_extreme != previous_extreme)
            {
                // if this is not the first shift and phenomenon is smooth enough
                if (previous_shift_extreme != NONE && (current_extreme_index - previous_shift_extreme_index) / dx >= smoothness)
                {
                    // we are in the presence of a PEAK
                    if (previous_shift_extreme == MIN && current_extreme == MIN)
                    {
                        // computing peak height
                        int peak_max_idx;
                        double peak_height = compute_peak_height(histogram, previous_shift_extreme_index, current_extreme_index, &peak_max_idx);

                        // adding peak statistics only if it has the required height
                        if (peak_height >= height_thres)
                        {
                            // computing peak area
                            double peak_area = compute_peak_area(histogram, previous_shift_extreme_index, current_extreme_index);

                            PeakInfo pi = PeakInfo(previous_shift_extreme_index, current_extreme_index, peak_max_idx, peak_height, peak_area);
                            peaks.push_back(pi);

                        }
                    }

                    // otherwise, it's not relevant to our case
                }

                // update previous shift extreme
                previous_shift_extreme = previous_extreme;
                previous_shift_extreme_index = previous_extreme_index;
            }

            // updating previous extreme point
            previous_extreme = current_extreme;
            previous_extreme_index = current_extreme_index;
        }

        // updating intervals
        b = c;
        c += offset;

        // updating last interval monotony
        ip_monotony = ic_monotony;
    }

    // returning the result
    return peaks;
}

/// <summary>
/// Classifies an extreme point according to monotony
/// </summary>
/// <param name="previous">Monotony of the previous interval</param>
/// <param name="current">Monotony of the current interval</param>
/// <returns>Extreme point type</returns>
ExtremeType get_extreme_type(Monotony previous, Monotony current)
{
    if (previous == GROW && current == STABLE ||
        previous == GROW && current == ABATE ||
        previous == STABLE && current == ABATE)
        return MAX;

    if (previous == ABATE && current == STABLE ||
        previous == ABATE && current == GROW ||
        previous == STABLE && current == GROW)
        return MIN;

    return NONE;
}

/// <summary>
/// Finds the monotony in an interval using central differences
/// </summary>
/// <param name="fa">Function value at a [f(a)]</param>
/// <param name="fb">Function value at b [f(b)]</param>
/// <param name="dx">Interval size</param>
/// <param name="grow_thres">Threshold for a growing function</param>
/// <param name="abate_thres">Threshold for an abating function</param>
/// <returns>The type of monotony</returns>
Monotony find_monotony(double fa, double fb, int dx, double growth_thres, double abate_thres)
{
    // computing the slope of the segment
    double m = (fb - fa) / dx;

    if (m >= growth_thres)
        return GROW;

    if (m <= abate_thres)
        return ABATE;

    return STABLE;
}

/// <summary>
/// Computes peak height
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <param name="peak_max_idx">Index at which peak reaches its max value</param>
/// <returns>Peak height</returns>
/// <remarks>Does not perform histogram normalization</remarks>
double compute_peak_height(std::vector<double> histogram, int lb, int ub, int *peak_max_idx)
{
    // getting extremes values
    double lb_value = histogram[lb];
    double ub_value = histogram[ub];

    // getting the lowest of both extremes
    double peak_base_value = std::min(lb_value, ub_value);

    // analyzing peak for extracting 'height'
    double peak_max_value = peak_base_value;
    int max_idx = lb;
    for (int i = lb; i <= ub; i++)
    {
        // updating peak max value and corresponding index
        if (histogram[i] > peak_max_value)
        {
            peak_max_value = histogram[i];
            max_idx = i;
        }
    }

    // storing peaks maximum intensity value index
    *peak_max_idx = max_idx;

    // computing peak 'height'
    return peak_max_value - peak_base_value;
}

/// <summary>
/// Computes peak area
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <returns>Peak area computed</returns>
/// <remarks>Does not perform histogram normalization</remarks>
double compute_peak_area(std::vector<double> histogram, int lb, int ub)
{
    // computing area of type A1 for now
    return compute_peak_area1(histogram, lb, ub);
}

// area of type A1
double compute_peak_area1(std::vector<double> histogram, int lb, int ub)
{
    // getting extremes values
    double lb_value = histogram[lb];
    double ub_value = histogram[ub];

    // getting the lowest of both extremes
    double peak_base_value = std::min(lb_value, ub_value);

    double area = 0.0;
    for (int i = lb; i <= ub; i++)
    {
        // computing area of type A1 (taking lowest extreme as baseline)
        area += (histogram[i] - peak_base_value);
    }

    // returning the area computed
    return area;
}

// area of type A2
double compute_peak_area2(std::vector<double> histogram, int lb, int ub)
{
    // getting extremes values
    double lb_value = histogram[lb];
    double ub_value = histogram[ub];

    // computing segment between peak extremes
    double m = (ub_value - lb_value) / (ub - lb);
    double n = ub_value - m * ub;   // evaluating in 'ub' without loss of generality

    double area = 0.0;
    for (int i = lb; i <= ub; i++)
    {
        // computing area of type A2 (taking segment between both extremes as baseline)
        area += histogram[i] - (m * i + n);
    }

    // returning the area computed
    return area;
}

/// <summary>
/// Computes peak statistics
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <returns>Peak statistics</returns>
PeakInfo compute_peak_statistics(std::vector<double> histogram, int lb, int ub)
{
    // computing peak height
    int peak_max_idx;
    double peak_height = compute_peak_height(histogram, lb, ub, &peak_max_idx);

    // computing peak area
    double peak_area = compute_peak_area(histogram, lb, ub);

    // creating a statistics structure (DEBUG purposes only)
    PeakInfo stats = PeakInfo(lb, ub, peak_max_idx, peak_height, peak_area);

    // returning the statistics computed for the given peak
    return stats;
}

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
API_FUNC(PeakInfo*) find_in_histogram(double* histogram, int h_length, int dx, int smoothness, double growth_angle, double abate_angle, double height_thres, int* peaks_count)
{
    // building the parameters for the function
    std::vector<double> h;
    for (int i = 0; i < h_length; i++)
        h.push_back(histogram[i]);

    // actually calling the function
    std::vector<PeakInfo> peaks = find_in_histogram(h, dx, smoothness, growth_angle, abate_angle, height_thres);

    // setting the amount of peaks found
    *peaks_count = (int)peaks.size();

    // building the result data
    PeakInfo* result = new PeakInfo[peaks.size()];

    // data copy from vector to pointer
    for (std::size_t i = 0; i < peaks.size(); i++)
        result[i] = peaks[i];

    // returning the list of peaks and corresponding info
    return result;
}

/// <summary>
/// Computes peak statistics for further analysis
/// </summary>
/// <param name="histogram">Histogram involved</param>
/// <param name="lb">Peak start index</param>
/// <param name="ub">Peak end index</param>
/// <returns>Peak height, area, etc.</returns>
/// <remarks>Does not perform histogram normalization</remarks>
API_FUNC(PeakInfo) compute_peak_statistics(double* histogram, int lb, int ub)
{
    // transforming the section of the histogram to a vector
    std::vector<double> h(histogram + lb, histogram + ub + 1);

    // computing peak statistics
    PeakInfo stats = compute_peak_statistics(h, 0, (int)h.size() - 1);

    // updating indexes (lb is the starting point)
    stats.lower_bound += lb;
    stats.upper_bound += lb;
    stats.height_index += lb;

    // returning computed peak statistics
    return stats;
}

/// <summary>
/// Frees a 'dangling' pointer of PeakInfo
/// </summary>
/// <param name="ptr">The pointer in question</param>
/// <remarks>This should be called after calling C-API 'find_in_histogram' function</remarks>
// deletes an already allocated (PeakInfo*)
API_FUNC(void) delete_peak_info_ptr(PeakInfo* ptr)
{
    delete[] ptr;
}
