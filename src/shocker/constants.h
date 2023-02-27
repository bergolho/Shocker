#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <vector>
#include <string>
#include <list>
#include <algorithm>

// ============================================================================================================================================
// CONSTANTS AND MACROS 
// ============================================================================================================================================
#define UM_TO_CM 0.0001
#define UM_TO_MM 0.001
#define CM_TO_UM 1.0E+05
#define MM_TO_UM 1.0E+03
#define M_TO_UM 1.0E+06
#define M_TO_MM 1.0E+03
#define MS_TO_S 0.001
#define S_TO_MS 1000.0
#define CM_S_TO_M_S 0.01
#define M_S_TO_UM_MS 1000.0
#define UM_MS_TO_M_S 0.001;
#define MS_TO_US 1000.0
#define CM_TO_MM 10.0
#define CM_TO_M 0.01
#define USE_ADJUST_DIAMETER false
#define USE_REGION_RADIUS false

#define GRAY 0
#define GREEN 1
#define RED 2
#define YELLOW 3
#define PURPLE 4
#define CYAN 5

static const uint32_t NTOSS = 10;                   // Number of tosses for a new terminal
//static const uint32_t NCONN = 400;                  // Number of segments to test for connection
static const double FACTOR = 0.95;                  // Reduction factor for the distance criterion
static const double L_MIN_LIMIT = 1.0e-08;          // Limit for the minimum distance value
// ============================================================================================================================================

#endif
