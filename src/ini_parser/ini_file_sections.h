//
// Created by sachetto on 12/10/17.
//

#ifndef MONOALG3D_INI_FILE_HEADERS_H
#define MONOALG3D_INI_FILE_HEADERS_H

#define MAIN_SECTION "main"
#define SAVE_NETWORK_SECTION "save_network"
#define CLOUD_SECTION "cloud_points"
#define LOCAL_OPT_SECTION "local_optimization"
#define PMJ_SECTION "pmj"
#define COST_FUNCTION_SECTION "cost_function"
#define PRUNING_SECTION "pruning"

#define MATCH_SECTION_AND_NAME(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
#define MATCH_SECTION(s) strcmp(section, s) == 0
#define MATCH_NAME(v) strcmp(name, v) == 0
#define SECTION_STARTS_WITH(s) strncmp(section, s, strlen(s)) == 0
#define A_STARTS_WITH_B(a ,b) strncmp(a, b, strlen(b)) == 0


#endif //MONOALG3D_INI_FILE_HEADERS_H
