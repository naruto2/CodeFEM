#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

__kernel void hello(__global char* string)
{
   string[0] = 'H';
   string[1] = 'e';
   string[2] = 'l';
   string[3] = 'l';
   string[4] = 'o';
   string[5] = ',';
   string[6] = ' ';
   string[7] = 'W';
   string[8] = 'o';
   string[9] = 'r';
   string[10] = 'l';
   string[11] = 'A'+get_global_size(0);
   string[12] = 'A'+get_global_id(0);
   string[13] = '\0';
}
