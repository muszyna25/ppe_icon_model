
#line 1 "mtime_iso8601.rl"
/*! \cond PRIVATE */
/**
 * @brief ISO 8601_2004 complaint Time.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note USAGE: Compile this rl file and generate iso8601.c file. 
	  Compile the iso8601.c file using gcc. iso8601.h file needs to be edited seperately. 
          match_found = 1 => DATE/DATETIME. match_found = 2 => Duration. Else non-compliant string and hence REJECT.
	  Due to application requirements, current implementation allows year in the range 2147483647 and -2147483648 only!
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <stdbool.h>

#include "mtime_iso8601.h"

#define MAX_BUFFER_LENGTH 132

#define YEAR_UPPER_BOUND 2147483647
#define YEAR_LOWER_BOUND -2147483648

/* Allowed year range = 2147483647 TO -2147483648  */
bool RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = false;




#line 38 "mtime_iso8601.c"
static const char _date_machine_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	5, 1, 6, 1, 7, 1, 8, 1, 
	9, 1, 10, 1, 11, 1, 12, 1, 
	13, 1, 14, 1, 15, 1, 16, 1, 
	17, 1, 18, 1, 19, 2, 0, 5, 
	2, 0, 8, 2, 0, 9, 2, 0, 
	10, 2, 0, 11, 2, 0, 12, 2, 
	6, 9, 2, 7, 10, 2, 8, 0, 
	2, 9, 0, 2, 9, 11, 2, 10, 
	0, 2, 10, 11, 2, 11, 0, 2, 
	11, 9, 2, 11, 10, 3, 0, 9, 
	10, 3, 0, 10, 11, 3, 0, 11, 
	9, 3, 0, 11, 10, 3, 3, 0, 
	4, 3, 7, 0, 10, 3, 7, 9, 
	10, 3, 7, 10, 9, 3, 9, 0, 
	11, 3, 9, 10, 11, 3, 9, 11, 
	10, 3, 10, 0, 11, 3, 10, 11, 
	0, 3, 10, 11, 9, 4, 0, 10, 
	11, 9, 4, 7, 0, 9, 10, 4, 
	9, 0, 10, 11
};

static const short _date_machine_key_offsets[] = {
	0, 0, 8, 12, 21, 24, 26, 30, 
	32, 34, 38, 40, 42, 47, 49, 54, 
	59, 62, 63, 66, 70, 71, 72, 76, 
	77, 78, 83, 84, 88, 92, 95, 101, 
	109, 123, 127, 129, 131, 133, 143, 149, 
	156, 166, 170, 176, 181, 187, 189, 192, 
	194, 197, 200, 203, 210, 214, 217, 219, 
	226, 230, 232, 236, 246, 249, 257, 260, 
	268, 271, 282, 285, 287, 292, 294, 299, 
	302, 306, 307, 311, 319, 327, 335, 343, 
	353, 360, 369, 376, 380, 388, 399, 402, 
	408, 415, 419, 431, 442, 450, 458, 468, 
	476, 484, 494, 509, 517, 527, 535, 545, 
	555, 558, 568, 579, 588, 597, 606, 622, 
	629, 640, 647, 662, 665, 671, 682, 686, 
	695, 706, 719, 734, 738, 749, 753, 763, 
	773, 776, 782, 791, 806, 817, 832, 840, 
	850, 858, 868, 877, 886, 891, 900, 906, 
	914, 925, 932, 941, 949, 952, 958, 963, 
	971, 974, 979, 982, 988, 990, 992, 994, 
	999, 1001, 1005, 1007, 1010, 1013, 1014, 1015, 
	1018, 1020, 1023, 1026, 1028, 1030, 1034, 1039, 
	1041, 1042, 1043, 1047, 1052, 1056, 1059, 1063, 
	1074, 1078, 1088, 1091, 1094, 1098, 1100, 1102, 
	1104, 1105, 1110, 1112, 1116, 1117, 1119, 1121, 
	1123, 1126, 1128, 1130, 1132, 1134, 1138, 1141, 
	1144, 1146, 1147, 1148, 1152, 1163, 1166, 1173, 
	1179, 1186, 1203, 1214, 1221, 1231, 1244, 1250, 
	1266, 1274, 1280, 1288, 1292, 1299, 1302, 1307, 
	1315, 1321, 1324, 1329, 1334, 1338, 1341, 1345, 
	1351, 1351, 1354
};

static const char _date_machine_trans_keys[] = {
	10, 36, 43, 45, 80, 84, 48, 57, 
	45, 80, 48, 57, 10, 32, 48, 49, 
	50, 51, 84, 9, 13, 50, 48, 49, 
	48, 57, 10, 58, 90, 122, 48, 53, 
	48, 57, 10, 58, 90, 122, 48, 53, 
	48, 57, 10, 44, 46, 90, 122, 48, 
	57, 10, 90, 122, 48, 57, 10, 90, 
	122, 48, 57, 10, 90, 122, 10, 52, 
	48, 51, 10, 58, 90, 122, 48, 48, 
	10, 58, 90, 122, 48, 48, 10, 44, 
	46, 90, 122, 48, 10, 48, 90, 122, 
	10, 48, 90, 122, 48, 49, 57, 10, 
	58, 90, 122, 49, 57, 10, 32, 50, 
	84, 9, 13, 48, 49, 10, 32, 45, 
	50, 58, 84, 90, 122, 9, 13, 48, 
	49, 51, 57, 48, 51, 49, 50, 49, 
	57, 48, 57, 48, 49, 10, 32, 50, 
	84, 9, 13, 48, 49, 51, 57, 10, 
	58, 90, 122, 48, 57, 10, 52, 58, 
	90, 122, 48, 51, 10, 32, 50, 51, 
	52, 84, 9, 13, 48, 49, 48, 50, 
	51, 57, 10, 58, 90, 122, 48, 57, 
	52, 48, 51, 53, 57, 10, 58, 90, 
	122, 48, 57, 48, 57, 54, 48, 53, 
	48, 54, 45, 48, 57, 45, 48, 57, 
	45, 48, 57, 45, 48, 49, 50, 51, 
	52, 57, 48, 49, 50, 51, 48, 49, 
	57, 49, 57, 10, 45, 50, 48, 49, 
	51, 57, 48, 50, 51, 57, 48, 57, 
	45, 48, 49, 57, 10, 45, 48, 58, 
	90, 122, 49, 53, 54, 57, 45, 48, 
	57, 10, 45, 90, 122, 48, 53, 54, 
	57, 45, 48, 57, 10, 45, 90, 122, 
	44, 46, 48, 57, 45, 48, 57, 10, 
	32, 45, 50, 84, 9, 13, 48, 49, 
	51, 57, 50, 48, 49, 48, 57, 10, 
	90, 122, 48, 53, 48, 57, 10, 90, 
	122, 48, 53, 52, 48, 51, 10, 48, 
	90, 122, 48, 10, 48, 90, 122, 10, 
	45, 90, 122, 48, 53, 54, 57, 10, 
	45, 90, 122, 48, 53, 54, 57, 10, 
	45, 90, 122, 44, 46, 48, 57, 10, 
	45, 90, 122, 48, 53, 54, 57, 10, 
	45, 52, 53, 90, 122, 48, 51, 54, 
	57, 10, 45, 48, 90, 122, 49, 57, 
	10, 45, 48, 90, 122, 44, 46, 49, 
	57, 10, 45, 48, 90, 122, 49, 57, 
	45, 48, 49, 57, 10, 45, 90, 122, 
	44, 46, 48, 57, 10, 32, 45, 50, 
	84, 9, 13, 48, 49, 51, 57, 45, 
	48, 57, 45, 52, 48, 51, 53, 57, 
	10, 45, 48, 90, 122, 49, 57, 45, 
	48, 49, 57, 45, 48, 49, 50, 51, 
	58, 90, 122, 52, 53, 54, 57, 10, 
	32, 45, 48, 49, 50, 84, 9, 13, 
	51, 57, 10, 45, 90, 122, 48, 53, 
	54, 57, 10, 45, 90, 122, 48, 53, 
	54, 57, 10, 45, 90, 122, 44, 46, 
	48, 53, 54, 57, 10, 45, 90, 122, 
	48, 53, 54, 57, 10, 45, 90, 122, 
	44, 46, 48, 57, 10, 45, 90, 122, 
	44, 46, 48, 53, 54, 57, 10, 32, 
	45, 50, 84, 90, 122, 9, 13, 48, 
	49, 51, 53, 54, 57, 10, 45, 90, 
	122, 48, 53, 54, 57, 10, 45, 90, 
	122, 44, 46, 48, 53, 54, 57, 10, 
	45, 90, 122, 48, 53, 54, 57, 10, 
	45, 90, 122, 44, 46, 48, 53, 54, 
	57, 10, 45, 90, 122, 44, 46, 48, 
	53, 54, 57, 45, 48, 57, 10, 45, 
	52, 53, 90, 122, 48, 51, 54, 57, 
	10, 45, 48, 90, 122, 44, 46, 49, 
	53, 54, 57, 10, 45, 48, 90, 122, 
	49, 53, 54, 57, 10, 45, 48, 90, 
	122, 44, 46, 49, 57, 10, 45, 48, 
	90, 122, 44, 46, 49, 57, 10, 32, 
	45, 50, 51, 52, 53, 84, 90, 122, 
	9, 13, 48, 49, 54, 57, 10, 45, 
	48, 90, 122, 49, 57, 10, 45, 48, 
	90, 122, 44, 46, 49, 53, 54, 57, 
	10, 45, 48, 90, 122, 49, 57, 10, 
	32, 45, 50, 84, 90, 122, 9, 13, 
	48, 49, 51, 53, 54, 57, 45, 48, 
	57, 45, 52, 48, 51, 53, 57, 10, 
	45, 48, 90, 122, 44, 46, 49, 53, 
	54, 57, 45, 48, 49, 57, 10, 45, 
	48, 90, 122, 44, 46, 49, 57, 10, 
	32, 45, 50, 84, 9, 13, 48, 49, 
	51, 57, 10, 32, 45, 50, 51, 52, 
	84, 9, 13, 48, 49, 53, 57, 10, 
	32, 45, 48, 49, 50, 84, 90, 122, 
	9, 13, 51, 53, 54, 57, 45, 48, 
	49, 57, 10, 45, 48, 90, 122, 44, 
	46, 49, 53, 54, 57, 45, 48, 49, 
	57, 10, 45, 90, 122, 44, 46, 48, 
	53, 54, 57, 10, 45, 90, 122, 44, 
	46, 48, 53, 54, 57, 45, 48, 57, 
	45, 52, 48, 51, 53, 57, 10, 45, 
	48, 90, 122, 44, 46, 49, 57, 10, 
	32, 45, 50, 84, 90, 122, 9, 13, 
	48, 49, 51, 53, 54, 57, 10, 32, 
	45, 50, 84, 9, 13, 48, 49, 51, 
	57, 10, 32, 45, 50, 84, 90, 122, 
	9, 13, 48, 49, 51, 53, 54, 57, 
	10, 45, 90, 122, 48, 53, 54, 57, 
	10, 45, 90, 122, 44, 46, 48, 53, 
	54, 57, 10, 45, 90, 122, 48, 53, 
	54, 57, 10, 45, 52, 53, 90, 122, 
	48, 51, 54, 57, 10, 45, 48, 90, 
	122, 44, 46, 49, 57, 10, 45, 48, 
	90, 122, 49, 53, 54, 57, 45, 48, 
	50, 51, 57, 10, 45, 58, 90, 122, 
	48, 53, 54, 57, 45, 52, 48, 51, 
	53, 57, 10, 45, 48, 58, 90, 122, 
	49, 57, 10, 32, 45, 48, 49, 50, 
	84, 9, 13, 51, 57, 10, 45, 48, 
	90, 122, 49, 57, 10, 45, 48, 90, 
	122, 49, 53, 54, 57, 10, 45, 90, 
	122, 44, 46, 48, 57, 45, 48, 57, 
	45, 54, 48, 53, 55, 57, 45, 48, 
	54, 55, 57, 10, 48, 49, 50, 51, 
	84, 52, 57, 89, 48, 57, 68, 77, 
	89, 48, 57, 89, 48, 57, 10, 48, 
	49, 50, 51, 84, 48, 57, 68, 77, 
	10, 84, 50, 48, 49, 51, 53, 48, 
	57, 46, 72, 77, 83, 48, 57, 83, 
	48, 57, 83, 48, 57, 83, 10, 10, 
	48, 53, 48, 57, 46, 77, 83, 10, 
	48, 53, 48, 57, 46, 83, 48, 51, 
	52, 57, 10, 51, 84, 48, 50, 48, 
	57, 68, 48, 48, 49, 50, 57, 89, 
	48, 49, 50, 57, 68, 89, 48, 57, 
	89, 48, 57, 48, 89, 49, 57, 10, 
	32, 48, 49, 50, 51, 84, 9, 13, 
	52, 57, 48, 89, 49, 57, 10, 48, 
	58, 68, 77, 89, 90, 122, 49, 57, 
	89, 48, 57, 10, 45, 84, 48, 49, 
	50, 51, 48, 57, 68, 77, 10, 84, 
	45, 50, 48, 49, 51, 53, 48, 57, 
	46, 72, 77, 83, 45, 10, 45, 48, 
	53, 48, 57, 46, 77, 83, 10, 45, 
	48, 53, 48, 57, 46, 83, 48, 51, 
	52, 57, 10, 45, 84, 51, 48, 50, 
	48, 57, 68, 48, 48, 49, 50, 57, 
	10, 32, 50, 84, 89, 9, 13, 48, 
	49, 51, 57, 89, 48, 57, 10, 58, 
	89, 90, 122, 48, 57, 52, 89, 48, 
	51, 53, 57, 10, 58, 89, 90, 122, 
	48, 57, 10, 32, 45, 50, 58, 68, 
	77, 84, 89, 90, 122, 9, 13, 48, 
	49, 51, 57, 10, 32, 50, 84, 89, 
	9, 13, 48, 49, 51, 57, 10, 58, 
	89, 90, 122, 48, 57, 10, 52, 58, 
	89, 90, 122, 48, 51, 53, 57, 10, 
	32, 50, 51, 52, 84, 89, 9, 13, 
	48, 49, 53, 57, 50, 89, 48, 49, 
	51, 57, 10, 32, 45, 50, 58, 68, 
	84, 89, 90, 122, 9, 13, 48, 49, 
	51, 57, 10, 58, 68, 89, 90, 122, 
	48, 57, 52, 89, 48, 51, 53, 57, 
	10, 58, 68, 89, 90, 122, 48, 57, 
	68, 89, 48, 57, 48, 54, 89, 49, 
	53, 55, 57, 89, 48, 57, 89, 48, 
	54, 55, 57, 45, 48, 49, 50, 51, 
	84, 52, 57, 48, 49, 50, 51, 52, 
	57, 89, 48, 57, 68, 77, 89, 48, 
	57, 89, 48, 49, 50, 57, 68, 89, 
	48, 57, 89, 48, 57, 48, 89, 49, 
	57, 45, 50, 48, 49, 51, 53, 50, 
	48, 49, 50, 48, 49, 0
};

static const char _date_machine_single_lengths[] = {
	0, 6, 2, 7, 1, 0, 4, 0, 
	0, 4, 0, 0, 5, 0, 3, 3, 
	3, 1, 1, 4, 1, 1, 4, 1, 
	1, 5, 1, 4, 4, 1, 4, 4, 
	8, 2, 0, 0, 0, 4, 4, 5, 
	6, 0, 4, 1, 4, 0, 1, 0, 
	1, 1, 1, 5, 4, 1, 0, 3, 
	0, 0, 2, 6, 1, 4, 1, 4, 
	1, 5, 1, 0, 3, 0, 3, 1, 
	4, 1, 4, 4, 4, 4, 4, 6, 
	5, 5, 5, 2, 4, 5, 1, 2, 
	5, 2, 8, 7, 4, 4, 4, 4, 
	4, 4, 7, 4, 4, 4, 4, 4, 
	1, 6, 5, 5, 5, 5, 10, 5, 
	5, 5, 7, 1, 2, 5, 2, 5, 
	5, 7, 9, 2, 5, 2, 4, 4, 
	1, 2, 5, 7, 5, 7, 4, 4, 
	4, 6, 5, 5, 1, 5, 2, 6, 
	7, 5, 5, 4, 1, 2, 1, 6, 
	1, 3, 1, 6, 0, 2, 2, 1, 
	0, 4, 0, 1, 1, 1, 1, 1, 
	0, 3, 1, 0, 2, 0, 3, 0, 
	1, 1, 0, 1, 2, 1, 2, 7, 
	2, 8, 1, 3, 4, 0, 2, 2, 
	1, 1, 0, 4, 1, 2, 0, 0, 
	3, 2, 0, 0, 2, 0, 3, 1, 
	0, 1, 1, 0, 5, 1, 5, 2, 
	5, 11, 5, 5, 6, 7, 2, 10, 
	6, 2, 6, 2, 3, 1, 1, 6, 
	4, 1, 3, 1, 2, 1, 2, 2, 
	0, 1, 1
};

static const char _date_machine_range_lengths[] = {
	0, 1, 1, 1, 1, 1, 0, 1, 
	1, 0, 1, 1, 0, 1, 1, 1, 
	0, 0, 1, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 1, 1, 2, 
	3, 1, 1, 1, 1, 3, 1, 1, 
	2, 2, 1, 2, 1, 1, 1, 1, 
	1, 1, 1, 1, 0, 1, 1, 2, 
	2, 1, 1, 2, 1, 2, 1, 2, 
	1, 3, 1, 1, 1, 1, 1, 1, 
	0, 0, 0, 2, 2, 2, 2, 2, 
	1, 2, 1, 1, 2, 3, 1, 2, 
	1, 1, 2, 2, 2, 2, 3, 2, 
	2, 3, 4, 2, 3, 2, 3, 3, 
	1, 2, 3, 2, 2, 2, 3, 1, 
	3, 1, 4, 1, 2, 3, 1, 2, 
	3, 3, 3, 1, 3, 1, 3, 3, 
	1, 2, 2, 4, 3, 4, 2, 3, 
	2, 2, 2, 2, 2, 2, 2, 1, 
	2, 1, 2, 2, 1, 2, 2, 1, 
	1, 1, 1, 0, 1, 0, 0, 2, 
	1, 0, 1, 1, 1, 0, 0, 1, 
	1, 0, 1, 1, 0, 2, 1, 1, 
	0, 0, 2, 2, 1, 1, 1, 2, 
	1, 1, 1, 0, 0, 1, 0, 0, 
	0, 2, 1, 0, 0, 0, 1, 1, 
	0, 0, 1, 1, 0, 2, 0, 1, 
	1, 0, 0, 2, 3, 1, 1, 2, 
	1, 3, 3, 1, 2, 3, 2, 3, 
	1, 2, 1, 1, 2, 1, 2, 1, 
	1, 1, 1, 2, 1, 1, 1, 2, 
	0, 1, 1
};

static const short _date_machine_index_offsets[] = {
	0, 0, 8, 12, 21, 24, 26, 31, 
	33, 35, 40, 42, 44, 50, 52, 57, 
	62, 66, 68, 71, 76, 78, 80, 85, 
	87, 89, 95, 97, 102, 107, 110, 116, 
	123, 135, 139, 141, 143, 145, 153, 159, 
	166, 175, 178, 184, 188, 194, 196, 199, 
	201, 204, 207, 210, 217, 222, 225, 227, 
	233, 236, 238, 242, 251, 254, 261, 264, 
	271, 274, 283, 286, 288, 293, 295, 300, 
	303, 308, 310, 315, 322, 329, 336, 343, 
	352, 359, 367, 374, 378, 385, 394, 397, 
	402, 409, 413, 424, 434, 441, 448, 456, 
	463, 470, 478, 490, 497, 505, 512, 520, 
	528, 531, 540, 549, 557, 565, 573, 587, 
	594, 603, 610, 622, 625, 630, 639, 643, 
	651, 660, 671, 684, 688, 697, 701, 709, 
	717, 720, 725, 733, 745, 754, 766, 773, 
	781, 788, 797, 805, 813, 817, 825, 830, 
	838, 848, 855, 863, 870, 873, 878, 882, 
	890, 893, 898, 901, 908, 910, 913, 916, 
	920, 922, 927, 929, 932, 935, 937, 939, 
	942, 944, 948, 951, 953, 956, 959, 964, 
	966, 968, 970, 973, 977, 981, 984, 988, 
	998, 1002, 1012, 1015, 1019, 1024, 1026, 1029, 
	1032, 1034, 1038, 1040, 1045, 1047, 1050, 1052, 
	1054, 1058, 1061, 1063, 1065, 1068, 1071, 1075, 
	1078, 1080, 1082, 1084, 1087, 1096, 1099, 1106, 
	1111, 1118, 1133, 1142, 1149, 1158, 1169, 1174, 
	1188, 1196, 1201, 1209, 1213, 1219, 1222, 1226, 
	1234, 1240, 1243, 1248, 1252, 1256, 1259, 1263, 
	1268, 1269, 1272
};

static const unsigned char _date_machine_indicies[] = {
	0, 2, 2, 3, 5, 6, 4, 1, 
	7, 8, 4, 1, 10, 9, 11, 12, 
	13, 14, 9, 9, 1, 16, 15, 1, 
	17, 1, 18, 19, 20, 20, 1, 21, 
	1, 22, 1, 18, 23, 20, 20, 1, 
	24, 1, 25, 1, 18, 26, 26, 20, 
	20, 1, 27, 1, 18, 20, 20, 28, 
	1, 18, 20, 20, 29, 1, 18, 20, 
	20, 1, 18, 1, 30, 17, 1, 18, 
	31, 20, 20, 1, 32, 1, 33, 1, 
	18, 34, 20, 20, 1, 35, 1, 36, 
	1, 18, 37, 37, 20, 20, 1, 38, 
	1, 18, 39, 20, 20, 1, 18, 29, 
	20, 20, 1, 40, 41, 1, 18, 19, 
	20, 20, 42, 1, 10, 9, 16, 9, 
	9, 15, 1, 10, 9, 43, 45, 19, 
	9, 20, 20, 9, 44, 42, 1, 46, 
	48, 47, 1, 49, 1, 49, 1, 49, 
	1, 10, 9, 51, 9, 9, 50, 17, 
	1, 18, 19, 20, 20, 17, 1, 18, 
	30, 19, 20, 20, 17, 1, 10, 9, 
	51, 17, 30, 9, 9, 50, 1, 41, 
	52, 1, 18, 19, 20, 20, 42, 1, 
	53, 52, 54, 1, 18, 31, 20, 20, 
	42, 1, 42, 1, 55, 54, 1, 42, 
	1, 56, 57, 1, 56, 58, 1, 56, 
	59, 1, 60, 61, 62, 63, 64, 65, 
	1, 66, 67, 68, 14, 1, 69, 70, 
	1, 42, 1, 18, 43, 45, 44, 42, 
	1, 70, 54, 1, 54, 1, 56, 71, 
	72, 1, 18, 56, 73, 19, 20, 20, 
	74, 75, 1, 56, 76, 1, 18, 56, 
	20, 20, 77, 65, 1, 56, 78, 1, 
	18, 56, 20, 20, 26, 65, 1, 56, 
	65, 1, 80, 79, 56, 82, 79, 79, 
	81, 76, 1, 84, 83, 1, 85, 1, 
	18, 20, 20, 86, 1, 87, 1, 18, 
	20, 20, 24, 1, 88, 85, 1, 18, 
	89, 20, 20, 1, 90, 1, 18, 35, 
	20, 20, 1, 18, 56, 20, 20, 91, 
	92, 1, 18, 56, 20, 20, 93, 78, 
	1, 18, 56, 20, 20, 26, 76, 1, 
	18, 56, 20, 20, 73, 65, 1, 18, 
	56, 94, 77, 20, 20, 91, 65, 1, 
	18, 56, 95, 20, 20, 78, 1, 18, 
	56, 96, 20, 20, 26, 65, 1, 18, 
	56, 97, 20, 20, 65, 1, 56, 98, 
	65, 1, 18, 56, 20, 20, 37, 65, 
	1, 80, 79, 56, 100, 79, 79, 99, 
	65, 1, 56, 92, 1, 56, 101, 92, 
	65, 1, 18, 56, 102, 20, 20, 65, 
	1, 56, 96, 65, 1, 56, 103, 104, 
	105, 106, 19, 20, 20, 74, 75, 1, 
	80, 79, 56, 107, 108, 109, 79, 79, 
	110, 1, 18, 56, 20, 20, 111, 92, 
	1, 18, 56, 20, 20, 112, 113, 1, 
	18, 56, 20, 20, 26, 114, 76, 1, 
	18, 56, 20, 20, 115, 78, 1, 18, 
	56, 20, 20, 26, 78, 1, 18, 56, 
	20, 20, 26, 77, 65, 1, 80, 79, 
	56, 117, 79, 20, 20, 79, 116, 111, 
	92, 1, 18, 56, 20, 20, 118, 119, 
	1, 18, 56, 20, 20, 26, 120, 76, 
	1, 18, 56, 20, 20, 121, 113, 1, 
	18, 56, 20, 20, 26, 115, 78, 1, 
	18, 56, 20, 20, 26, 122, 65, 1, 
	56, 113, 1, 18, 56, 123, 112, 20, 
	20, 118, 113, 1, 18, 56, 124, 20, 
	20, 26, 114, 76, 1, 18, 56, 125, 
	20, 20, 115, 78, 1, 18, 56, 126, 
	20, 20, 26, 78, 1, 18, 56, 98, 
	20, 20, 26, 65, 1, 80, 79, 56, 
	117, 111, 127, 122, 79, 20, 20, 79, 
	116, 65, 1, 18, 56, 128, 20, 20, 
	113, 1, 18, 56, 129, 20, 20, 26, 
	77, 65, 1, 18, 56, 126, 20, 20, 
	78, 1, 80, 79, 56, 131, 79, 20, 
	20, 79, 130, 122, 65, 1, 56, 119, 
	1, 56, 132, 119, 113, 1, 18, 56, 
	133, 20, 20, 26, 77, 65, 1, 56, 
	134, 78, 1, 18, 56, 97, 20, 20, 
	26, 65, 1, 80, 79, 56, 109, 79, 
	79, 108, 110, 1, 80, 79, 56, 109, 
	110, 135, 79, 79, 108, 136, 1, 80, 
	79, 56, 137, 138, 139, 79, 20, 20, 
	79, 77, 65, 1, 56, 140, 141, 1, 
	18, 56, 142, 20, 20, 26, 73, 65, 
	1, 56, 143, 76, 1, 18, 56, 20, 
	20, 37, 77, 65, 1, 18, 56, 20, 
	20, 26, 73, 65, 1, 56, 141, 1, 
	56, 144, 141, 78, 1, 18, 56, 102, 
	20, 20, 26, 65, 1, 80, 79, 56, 
	139, 79, 20, 20, 79, 138, 77, 65, 
	1, 80, 79, 56, 82, 79, 79, 145, 
	76, 1, 80, 79, 56, 147, 79, 20, 
	20, 79, 146, 91, 92, 1, 18, 56, 
	20, 20, 148, 141, 1, 18, 56, 20, 
	20, 26, 149, 76, 1, 18, 56, 20, 
	20, 114, 76, 1, 18, 56, 150, 93, 
	20, 20, 148, 78, 1, 18, 56, 151, 
	20, 20, 26, 76, 1, 18, 56, 129, 
	20, 20, 77, 65, 1, 56, 72, 152, 
	1, 18, 56, 19, 20, 20, 74, 75, 
	1, 56, 153, 152, 154, 1, 18, 56, 
	155, 31, 20, 20, 75, 1, 80, 79, 
	56, 156, 99, 100, 79, 79, 65, 1, 
	18, 56, 157, 20, 20, 92, 1, 18, 
	56, 158, 20, 20, 73, 65, 1, 18, 
	56, 20, 20, 37, 76, 1, 56, 75, 
	1, 56, 159, 154, 65, 1, 56, 75, 
	65, 1, 0, 160, 161, 162, 163, 165, 
	164, 1, 167, 166, 1, 169, 170, 167, 
	168, 1, 167, 168, 1, 0, 171, 172, 
	173, 174, 165, 1, 175, 1, 169, 170, 
	1, 0, 165, 1, 177, 176, 178, 1, 
	179, 1, 180, 181, 182, 183, 1, 184, 
	1, 183, 185, 1, 183, 186, 1, 183, 
	1, 0, 1, 0, 178, 1, 187, 1, 
	180, 182, 183, 1, 0, 188, 1, 189, 
	1, 180, 183, 1, 179, 187, 1, 0, 
	174, 165, 173, 1, 190, 1, 169, 1, 
	190, 1, 175, 190, 1, 167, 166, 191, 
	1, 169, 167, 168, 1, 167, 191, 1, 
	191, 167, 168, 1, 10, 9, 192, 193, 
	194, 195, 9, 9, 196, 1, 197, 199, 
	198, 1, 18, 200, 19, 202, 203, 199, 
	20, 20, 201, 1, 199, 200, 1, 0, 
	204, 6, 1, 205, 206, 207, 208, 1, 
	209, 1, 202, 203, 1, 0, 6, 1, 
	210, 1, 212, 211, 213, 1, 214, 1, 
	215, 216, 217, 183, 1, 180, 1, 0, 
	218, 1, 213, 1, 219, 1, 215, 217, 
	183, 1, 0, 220, 1, 221, 1, 222, 
	1, 215, 183, 1, 214, 219, 1, 0, 
	223, 6, 1, 208, 207, 1, 224, 1, 
	202, 1, 224, 1, 209, 224, 1, 10, 
	9, 226, 9, 199, 9, 225, 200, 1, 
	199, 227, 1, 18, 19, 199, 20, 20, 
	200, 1, 228, 199, 227, 200, 1, 18, 
	31, 199, 20, 20, 200, 1, 10, 9, 
	43, 230, 19, 202, 203, 9, 199, 20, 
	20, 9, 229, 201, 1, 10, 9, 232, 
	9, 199, 9, 231, 227, 1, 18, 19, 
	199, 20, 20, 227, 1, 18, 228, 19, 
	199, 20, 20, 227, 200, 1, 10, 9, 
	232, 227, 228, 9, 199, 9, 231, 200, 
	1, 233, 199, 198, 234, 1, 10, 9, 
	43, 230, 19, 202, 9, 199, 20, 20, 
	9, 229, 201, 1, 18, 19, 202, 199, 
	20, 20, 201, 1, 235, 199, 234, 236, 
	1, 18, 31, 202, 199, 20, 20, 201, 
	1, 202, 199, 201, 1, 236, 238, 199, 
	237, 200, 1, 199, 201, 1, 199, 201, 
	200, 1, 239, 160, 161, 162, 163, 240, 
	164, 1, 241, 242, 243, 244, 196, 1, 
	199, 245, 1, 202, 203, 199, 200, 1, 
	199, 245, 246, 1, 202, 199, 200, 1, 
	199, 246, 1, 246, 199, 200, 1, 210, 
	177, 176, 178, 1, 1, 16, 15, 1, 
	84, 83, 1, 0
};

static const unsigned char _date_machine_trans_targs[] = {
	240, 0, 2, 183, 48, 231, 192, 3, 
	151, 4, 241, 29, 41, 43, 46, 5, 
	18, 6, 240, 7, 17, 8, 9, 10, 
	11, 12, 13, 14, 15, 16, 19, 20, 
	21, 22, 23, 24, 25, 26, 27, 28, 
	30, 32, 31, 33, 37, 40, 34, 35, 
	36, 31, 38, 39, 42, 44, 45, 47, 
	3, 49, 50, 51, 52, 58, 140, 142, 
	149, 64, 53, 56, 57, 54, 55, 59, 
	90, 60, 65, 85, 61, 62, 63, 66, 
	242, 75, 79, 67, 71, 68, 69, 70, 
	72, 73, 74, 76, 78, 77, 80, 81, 
	82, 83, 84, 86, 87, 88, 89, 91, 
	120, 121, 132, 92, 98, 110, 114, 93, 
	94, 97, 95, 96, 99, 105, 100, 103, 
	101, 102, 104, 106, 107, 108, 109, 111, 
	112, 113, 115, 116, 117, 118, 119, 122, 
	131, 123, 128, 129, 124, 127, 125, 126, 
	130, 133, 134, 137, 135, 136, 138, 139, 
	141, 143, 148, 144, 145, 146, 147, 150, 
	152, 179, 181, 182, 154, 159, 153, 155, 
	154, 158, 174, 156, 178, 175, 177, 157, 
	160, 173, 168, 161, 162, 167, 170, 166, 
	163, 164, 165, 169, 171, 172, 176, 180, 
	184, 222, 225, 228, 186, 185, 217, 187, 
	186, 212, 191, 206, 188, 189, 211, 208, 
	210, 190, 193, 194, 205, 199, 195, 196, 
	197, 201, 198, 200, 202, 203, 204, 207, 
	209, 213, 215, 214, 216, 218, 221, 219, 
	220, 223, 224, 226, 227, 229, 230, 232, 
	239, 233, 235, 237, 238, 234, 236
};

static const unsigned char _date_machine_trans_actions[] = {
	5, 0, 101, 37, 1, 0, 0, 37, 
	0, 0, 3, 1, 1, 1, 1, 1, 
	1, 15, 3, 0, 23, 1, 17, 0, 
	1, 19, 0, 52, 21, 21, 15, 0, 
	1, 17, 0, 1, 19, 0, 52, 21, 
	15, 55, 13, 0, 61, 61, 1, 1, 
	1, 11, 43, 43, 15, 15, 0, 0, 
	7, 0, 0, 7, 0, 1, 1, 1, 
	1, 0, 1, 1, 1, 0, 9, 15, 
	55, 1, 61, 13, 17, 1, 19, 0, 
	3, 46, 46, 1, 1, 15, 1, 17, 
	15, 1, 17, 64, 15, 49, 43, 76, 
	17, 1, 19, 1, 1, 15, 1, 40, 
	40, 40, 40, 85, 146, 146, 109, 64, 
	89, 73, 70, 49, 43, 43, 151, 121, 
	46, 129, 1, 141, 70, 97, 76, 43, 
	133, 46, 1, 1, 137, 1, 82, 113, 
	58, 1, 1, 1, 125, 67, 1, 73, 
	79, 105, 43, 43, 117, 46, 93, 70, 
	15, 15, 0, 61, 46, 64, 49, 0, 
	1, 1, 1, 1, 1, 0, 0, 25, 
	0, 29, 27, 1, 1, 1, 1, 0, 
	1, 1, 1, 0, 0, 31, 33, 35, 
	0, 0, 0, 0, 1, 0, 0, 0, 
	1, 1, 1, 1, 1, 15, 55, 25, 
	0, 13, 29, 27, 0, 1, 1, 1, 
	1, 0, 0, 1, 1, 1, 0, 0, 
	31, 33, 0, 0, 0, 1, 0, 0, 
	0, 1, 1, 15, 15, 61, 61, 43, 
	43, 55, 15, 15, 0, 0, 0, 0, 
	0, 1, 1, 1, 1, 0, 0
};

static const int date_machine_start = 1;
static const int date_machine_first_final = 240;
static const int date_machine_error = 0;

static const int date_machine_en_main = 1;


#line 37 "mtime_iso8601.rl"



struct internal_datetime
  {
    char            sign_of_year;
    int64_t         year;
    int             month;
    int             day;
    char            flag_365type_day;
    int             hour;
    int             minute;
    int             second;
    int 	    ms;
    char	    timeZoneOffSet;
  };


static 
void 
date_machine( char *str, ISO8601_STATUS* stat, struct internal_datetime* dtObj, struct iso8601_duration* duObj)
  {
    char *p = str, *pe = str + strlen( str );
    char *ts, *te = 0;
    int cs;

    
#line 637 "mtime_iso8601.c"
	{
	cs = date_machine_start;
	}

#line 642 "mtime_iso8601.c"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_keys = _date_machine_trans_keys + _date_machine_key_offsets[cs];
	_trans = _date_machine_index_offsets[cs];

	_klen = _date_machine_single_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _date_machine_range_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _date_machine_indicies[_trans];
	cs = _date_machine_trans_targs[_trans];

	if ( _date_machine_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _date_machine_actions + _date_machine_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 65 "mtime_iso8601.rl"
	{
	    ts = p;
	  }
	break;
	case 1:
#line 70 "mtime_iso8601.rl"
	{
	    *stat = DATETIME_MATCH;
	  }
	break;
	case 2:
#line 75 "mtime_iso8601.rl"
	{
	    *stat = DURATION_MATCH;
	  }
	break;
	case 3:
#line 80 "mtime_iso8601.rl"
	{
	    dtObj->sign_of_year = (*p);
	  }
	break;
	case 4:
#line 85 "mtime_iso8601.rl"
	{
	    duObj->sign = (*p);
	  }
	break;
	case 5:
#line 90 "mtime_iso8601.rl"
	{
	    te = p+1; 
	    /* Reset ts to point to begining of string. */
	    ts = str;                      

	    char _year[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _year, ts, te-ts);
	    _year[MAX_BUFFER_LENGTH-1] = '\0';
	    
	    /* To ensure strtol works. */
	    if ( _year[0] == '$')
	      _year[0] = '-';

	    char *end;
	    dtObj->year = strtol(_year,&end, 10);
	    
	    /* Year might be too large. */
	    if ((errno == ERANGE) || (dtObj->year > YEAR_UPPER_BOUND) || (dtObj->year < YEAR_LOWER_BOUND))
	      {
		RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = true;
	      }
	  }
	break;
	case 6:
#line 114 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _month[3] = {0};
	    strncpy( _month, ts, te-ts);
	    _month[2] = '\0';
	    dtObj->month = atoi(_month);
	  }
	break;
	case 7:
#line 123 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _day[3] = {0};
	    strncpy( _day, ts, te-ts);
	    _day[2] = '\0';
	    dtObj->day = atoi(_day);
	    dtObj->flag_365type_day = 0;
	  }
	break;
	case 8:
#line 133 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _day[4] = {0};
	    strncpy( _day, ts, te-ts);
	    _day[3] = '\0';
	    dtObj->day = atoi(_day);
	    dtObj->flag_365type_day = 1;
	  }
	break;
	case 9:
#line 143 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _hour[3] = {0};
	    strncpy( _hour, ts, te-ts);
	    _hour[2] = '\0';
	    dtObj->hour = atoi(_hour);
	  }
	break;
	case 10:
#line 152 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _minute[3] = {0};
	    strncpy( _minute, ts, te-ts);
	    _minute[2] = '\0';
	    dtObj->minute = atoi(_minute);
	  }
	break;
	case 11:
#line 161 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _second[3] = {0};
	    strncpy( _second, ts, te-ts);
	    _second[2] = '\0';
	    dtObj->second = atoi(_second);                
	  }
	break;
	case 12:
#line 170 "mtime_iso8601.rl"
	{
	    te = p+1;
	    char _ms[4] = {0};
	    strncpy( _ms, ts, te-ts);
	    _ms[3] = '\0';
	    dtObj->ms = atoi(_ms);
	  }
	break;
	case 13:
#line 179 "mtime_iso8601.rl"
	{
	    dtObj->timeZoneOffSet = (*p);
	  }
	break;
	case 14:
#line 184 "mtime_iso8601.rl"
	{
	    te = p;                 
	    char _du_year[MAX_BUFFER_LENGTH] = {'\0'};
	    strncpy( _du_year, ts, te-ts);
	    _du_year[MAX_BUFFER_LENGTH-1] = '\0';

	    char *end;
            duObj->year = strtol(_du_year,&end, 10);

            if ((errno == ERANGE) || ((duObj->sign == '$') && (duObj->year > (-1*YEAR_LOWER_BOUND))) || ((duObj->year > YEAR_UPPER_BOUND) && (duObj->sign != '$')))
                RAISE_YEAR_OUT_OF_BOUND_EXCEPTION = true;

	  }
	break;
	case 15:
#line 199 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_month[3] = {0};
	    strncpy( _du_month, ts, te-ts);
	    _du_month[2] = '\0';
	    duObj->month = atoi(_du_month);
	  }
	break;
	case 16:
#line 208 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_day[3] = {0};
	    strncpy( _du_day, ts, te-ts);
	    _du_day[2] = '\0';
	    duObj->day = atoi(_du_day);     
	  }
	break;
	case 17:
#line 217 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_hour[3] = {0};
	    strncpy( _du_hour, ts, te-ts);
	    _du_hour[2] = '\0';
	    duObj->hour = atoi(_du_hour);
	  }
	break;
	case 18:
#line 226 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_minute[3] = {0};
	    strncpy( _du_minute, ts, te-ts);
	    _du_minute[2] = '\0';
	    duObj->minute = atoi(_du_minute);
	  }
	break;
	case 19:
#line 235 "mtime_iso8601.rl"
	{
	    te = p;
	    char _du_second[13] = {0};
	    char* _du_ms;
	    int _ms;
	    
	    strncpy( _du_second, ts, te-ts);
	    _du_second[12] = '\0';
	    duObj->second = atoi(_du_second);
	    
	    if(strstr(_du_second,"."))
	      {
		_du_ms = (strstr(_du_second,".")+1);
		if(_du_ms[0] == '-')
		  _du_ms = _du_ms + 1;

		_ms = atoi(_du_ms);

		if(strlen(_du_ms) == 1)
		  duObj->ms = _ms*100;
		else if(strlen(_du_ms) == 2)
		  duObj->ms = _ms*10;
		else
		  duObj->ms = _ms;        
	      }               
	  }
	break;
#line 934 "mtime_iso8601.c"
		}
	}

_again:
	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	_out: {}
	}

#line 549 "mtime_iso8601.rl"

    
}

/* Internal function which calls the state machine. */
static 
ISO8601_STATUS 
get_date_time(const char* buffer, struct iso8601_datetime* datetimeObj, struct iso8601_duration* durationObj)
  {
    /* Create a local buffer and copy the string to be tested. */
    char buf[MAX_BUFFER_LENGTH] = {'\0'};
    strncpy(buf,buffer,MAX_BUFFER_LENGTH);
    buf[MAX_BUFFER_LENGTH-1] = '\0';

    /* Initialize Success or Failure flag. */
    ISO8601_STATUS stat = FAILURE;

    /* Placeholder for values of DateTime and Duration. */
    struct internal_datetime dtObj = {0};
    struct iso8601_duration duObj = {0};

    /* Initialize month and day to 1. In case these values are not specified in the buffer
       the default value should be 1. For eg. Date 1999-12 should return Year 1999, Month
       12 and DATE as "1".
    */
    dtObj.month = 1;
    dtObj.day = 1;

    /* Ragel expects \n at the end of the string. */
    char* replace = strchr(buf, '\0');
    *replace = '\n';

    /*The fact that '-' sign is used to denote negative years as well as a seperator in datetime 
      causes the regex to fail in certain scenarios. The fix ( hack? ) is to replace the - sign 
      with a '$' sign and copy back the value after processing. 
    */
    if(buf[0] == '-')
      buf[0] = '$';


    /* Execute Ragel Machine. */
    date_machine(buf,&stat,&dtObj,&duObj);

   
    /* stat contains the type of match. */
    if((stat == DATETIME_MATCH) && (RAISE_YEAR_OUT_OF_BOUND_EXCEPTION == false))
      {       
	/* Set sign of year. */
	if(dtObj.sign_of_year == '$')
	  datetimeObj->sign_of_year = '-';
	else
	  datetimeObj->sign_of_year = '+';

	/* Year. */
	datetimeObj->year = dtObj.year;
	
	/* If this is a 365 type day, set month to 0. */
	if(dtObj.flag_365type_day == 1)
	  {
	    datetimeObj->month  = 0;
	  }
	else
	  {
	    datetimeObj->month  = dtObj.month;
	  }

	/* Set rest. */
	datetimeObj->day 	= dtObj.day;               
	datetimeObj->hour 	= dtObj.hour;
	datetimeObj->minute 	= dtObj.minute;
	datetimeObj->second 	= dtObj.second;
	datetimeObj->ms 	= dtObj.ms;
	datetimeObj->timeZoneOffSet = dtObj.timeZoneOffSet;
      }               
    else if((stat == DURATION_MATCH) && (RAISE_YEAR_OUT_OF_BOUND_EXCEPTION == false))
      {
	/* Set sign of duration. */
	if (duObj.sign == '$')
	  durationObj->sign = '-';
	else
	  durationObj->sign = '+';

	/* Set rest. */
	durationObj->year   = duObj.year;
	durationObj->month  = duObj.month;
	durationObj->day    = duObj.day;
	durationObj->hour   = duObj.hour;
	durationObj->minute = duObj.minute;
	durationObj->second = duObj.second;
	durationObj->ms     = duObj.ms;
      }
    else
      {
        stat = FAILURE;
      }

    return stat;
  }


/*Check DateTime string compliance and get DateTime values. */
ISO8601_STATUS 
verify_string_datetime(const char* test_string,struct iso8601_datetime* dummy_isoDtObj)
  {
    ISO8601_STATUS stat = FAILURE;
    struct iso8601_duration* dummy_isoDObj = new_iso8601_duration('+',0,0,0,0,0,0,0);
    if (dummy_isoDObj == NULL)
      return FAILURE;

    stat = get_date_time(test_string, dummy_isoDtObj, dummy_isoDObj);

    deallocate_iso8601_duration(dummy_isoDObj);

    return stat;
  }

/*Check TimeDelta string compliance and get duration values.*/
ISO8601_STATUS
verify_string_duration(const char* test_string, struct iso8601_duration* dummy_isoDObj)
  {
    ISO8601_STATUS stat = FAILURE;
    struct iso8601_datetime* dummy_isoDtObj = new_iso8601_datetime('+',0,0,0,0,0,0,0,'Z');
    if ( dummy_isoDtObj == NULL)
      return FAILURE;

    stat = get_date_time(test_string, dummy_isoDtObj, dummy_isoDObj);

    deallocate_iso8601_datetime(dummy_isoDtObj);

    return stat;
  }


struct iso8601_datetime* 
new_iso8601_datetime( char _sign_of_year, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms, char _timeZoneOffSet )
  {
    struct iso8601_datetime* isoDtObj = (struct iso8601_datetime*)calloc(1,sizeof(struct iso8601_datetime));
    if (isoDtObj == NULL)
      return NULL;

    isoDtObj->sign_of_year 	= _sign_of_year;
    isoDtObj->year 		= _year;
    isoDtObj->month 		= _month;
    isoDtObj->day 		= _day;
    isoDtObj->hour 		= _hour;
    isoDtObj->minute 		= _minute;
    isoDtObj->second 		= _second;
    isoDtObj->ms 		= _ms;
    isoDtObj->timeZoneOffSet 	= _timeZoneOffSet;

    return isoDtObj;
  }


void 
deallocate_iso8601_datetime(struct iso8601_datetime* iso8601_datetimeObj)
  {
    if ( iso8601_datetimeObj != NULL)
      {
	free(iso8601_datetimeObj);
	iso8601_datetimeObj = NULL;
      }
  }

struct iso8601_duration* 
new_iso8601_duration( char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms )
  {
    struct iso8601_duration* isoDObj = (struct iso8601_duration*)calloc(1,sizeof(struct iso8601_duration));
    if (isoDObj == NULL)
      return NULL;

    isoDObj->sign 	= _sign;
    isoDObj->year 	= _year;
    isoDObj->month 	= _month;
    isoDObj->day 	= _day;
    isoDObj->hour 	= _hour;
    isoDObj->minute 	= _minute;
    isoDObj->second 	= _second;
    isoDObj->ms 	= _ms;

    return isoDObj;
  }

void 
deallocate_iso8601_duration(struct iso8601_duration* iso8601_durationObj)
  {
    if ( iso8601_durationObj != NULL )
      {
	free(iso8601_durationObj);
	iso8601_durationObj = NULL;
      }
  }

/*! \endcond */
