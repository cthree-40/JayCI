/* iminmax.c */
/* 
 * Subfunctions to return minimum and maximum integers. 
 */

/* int_max: return maximum */
int int_max(int a, int b)
{
	int max;
	if (a >= b) {
		max = a;
		return max;
	} else {
		max = b;
		return max;
	}
}
/* int_min: return minimum */
int int_min(int a, int b)
{
	int min;
	if (a >= b) {
		min = b;
		return min;
	} else {
		min = a;
		return min;
	}
};
