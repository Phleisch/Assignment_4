__kernel void hello_world() {
	int index = get_local_id(0) + 4*get_local_id(1) +\
		4*4*get_local_id(2) + 4*4*4 * get_group_id(0);
	printf("Hello World! My threadId is %d\n", index);
}
