// Template file for the OpenCL Assignment 4

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>		// For use of gettimeofday function
#include <CL/cl.h>

// Define number of work items and work groups
#define NUM_WORK_ITEMS_PER_GROUP 256
//#define ARRAY_SIZE 10000

#define KERNEL_FILEPATH "./kernel.cl"

// This is a macro for checking the error variable.
#define CHK_ERROR(err) if ((err) != CL_SUCCESS) fprintf(stderr,"Error: %s\n",clGetErrorString(err));

// A errorCode to string converter (forward declaration)
const char* clGetErrorString(int);

cl_program compileKernelBoilerplate(
	char *kernelFilepath, cl_context context, cl_device_id *device_list) {
	cl_int err;
	char *sourceCode = 0;
	long fileLength;
	FILE *kernelFile = fopen (kernelFilepath, "rb");

	// File does not exist or is not accessible
	if (kernelFile == NULL) {
		fprintf(stderr, "Could not open %s\n", kernelFilepath);
		return NULL;
	}


	fseek (kernelFile, 0, SEEK_END);
	fileLength = ftell (kernelFile);
	fseek (kernelFile, 0, SEEK_SET);
	sourceCode = (char *) malloc(fileLength);

	// Could not allocate space for sourceCode
	if (sourceCode == NULL) {
		fprintf(stderr, "Malloc failed\n");
		return NULL;
	}

	fread (sourceCode, 1, fileLength, kernelFile);
	fclose (kernelFile);

	/* Create the OpenCL program */
	cl_program program = clCreateProgramWithSource(
			context, 1, (const char**) &sourceCode, NULL, &err);
	CHK_ERROR(err);

	err = clBuildProgram(program, 1, device_list, NULL, NULL, NULL);
	CHK_ERROR(err);

	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];
		clGetProgramBuildInfo(
				program, device_list[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer),
				buffer, &len);
		fprintf(stderr, "Build error: %s\n", buffer);
		exit(0);
	}

	return program;
}

/**
 * Return a timestamp with double precision.
 */
double cpuSecond() {
	struct timeval tp;
	gettimeofday(&tp,NULL);
	return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

int assert_equal(float *A, float *B, size_t array_size) {
  for (int i = 0; i < array_size; ++i) {
    if (abs(A[i] - B[i]) > 0.00001) {
      printf("Not equal in index %d: %f vs %f\n", i, A[i], B[i]);
      return 0;
    }
  }
  return 1;
}

void cpu_SAXPY(float a, float *X, float *Y, size_t array_size) {
  for (int i = 0; i < array_size; ++i) {
    Y[i] += a * X[i];
  }
}

int main(int argc, char **argv) {
#ifndef ARRAY_SIZE
  if (argc != 2)
    printf("Usage: %s <ARRAY_SIZE>\nOr define ARRAY_SIZE.\n", argv[0]);
  size_t ARRAY_SIZE = atoi(argv[1]);
  printf("%ld,", ARRAY_SIZE);
#endif
  cl_platform_id * platforms; cl_uint     n_platform;

  // Find OpenCL Platforms
  cl_int err = clGetPlatformIDs(0, NULL, &n_platform);
  CHK_ERROR(err);
  platforms = (cl_platform_id *) malloc(sizeof(cl_platform_id)*n_platform);
  err = clGetPlatformIDs(n_platform, platforms, NULL);
  CHK_ERROR(err);

  // Find and sort devices
  cl_device_id *device_list; cl_uint n_devices;
  err = clGetDeviceIDs( platforms[0], CL_DEVICE_TYPE_GPU, 0,NULL, &n_devices);
  CHK_ERROR(err);
  device_list = (cl_device_id *) malloc(sizeof(cl_device_id)*n_devices);
  err = clGetDeviceIDs( platforms[0],CL_DEVICE_TYPE_GPU, n_devices,
      device_list, NULL);
  CHK_ERROR(err);
  
  // Create and initialize an OpenCL context
  cl_context context =
    clCreateContext( NULL, n_devices, device_list, NULL, NULL, &err);
  CHK_ERROR(err);

  // Create a command queue
  cl_command_queue cmd_queue = clCreateCommandQueue(context, device_list[0],
      CL_QUEUE_PROFILING_ENABLE, &err);
  CHK_ERROR(err); 

  /* Create the OpenCL program */
  cl_program program = compileKernelBoilerplate(KERNEL_FILEPATH, context,
      device_list);

  /* Create a kernel object referencing our "SAXPY" kernel */
  cl_kernel kernel = clCreateKernel(program, "SAXPY", &err);
  CHK_ERROR(err);

  size_t n_workitem = ((ARRAY_SIZE + (NUM_WORK_ITEMS_PER_GROUP - 1)) / NUM_WORK_ITEMS_PER_GROUP) * NUM_WORK_ITEMS_PER_GROUP;
  size_t workgroup_size = NUM_WORK_ITEMS_PER_GROUP;

  /* Host arguments */
  size_t byte_size = ARRAY_SIZE * sizeof(float);
  float X[ARRAY_SIZE], Y[ARRAY_SIZE], Y_reference[ARRAY_SIZE];
  float a = 1.5f;
  for (int i = 0; i < ARRAY_SIZE; ++i) {
    X[i] = 2.f * i;
    Y[i] = -1.f * i;
    Y_reference[i] = -1.f * i;
  }

  /* Allocate device buffers */
  cl_mem X_dev =
    clCreateBuffer(context, CL_MEM_READ_ONLY, byte_size, NULL, &err);
  CHK_ERROR(err);
  cl_mem Y_dev =
    clCreateBuffer(context, CL_MEM_READ_WRITE, byte_size, NULL, &err);
  CHK_ERROR(err);

  /* Write device buffers */
  err = clEnqueueWriteBuffer(cmd_queue, X_dev, CL_TRUE, 0, byte_size, X, 0,
      NULL, NULL);
  CHK_ERROR(err);
  err = clEnqueueWriteBuffer(cmd_queue, Y_dev, CL_TRUE, 0, byte_size, Y, 0,
      NULL, NULL);
  CHK_ERROR(err);

  /* Compute CPU result */
	double startTime = cpuSecond();
#ifdef ARRAY_SIZE
  printf("Computing SAXPY on the CPU... ");
#endif
  cpu_SAXPY(a, X, Y_reference, ARRAY_SIZE);
#ifdef ARRAY_SIZE
	printf("Done in %f seconds.\n", cpuSecond() - startTime);
#else
	printf("%f,", cpuSecond() - startTime);
#endif

  /* Enqueue the kernel */
  cl_ulong array_size = ARRAY_SIZE;
  err = clSetKernelArg(kernel, 0, sizeof(float), (void *) &a);
  CHK_ERROR(err);
  err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &X_dev);
  CHK_ERROR(err);
  err = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &Y_dev);
  CHK_ERROR(err);
  err = clSetKernelArg(kernel, 3, sizeof(cl_ulong), (void *) &array_size);
  CHK_ERROR(err);

  cl_event event;
  err = clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL, &n_workitem, 
      &workgroup_size, 0, NULL, &event);
  CHK_ERROR(err);

  /* Get result back to the host */
  err = clEnqueueReadBuffer(cmd_queue, Y_dev, CL_TRUE, 0, byte_size, Y, 0,
      NULL, NULL);
  CHK_ERROR(err);

  /* Wait and make sure everything finished */
#ifdef ARRAY_SIZE
  printf("Computing SAXPY on the GPU... ");
#endif
  err = clFlush(cmd_queue);
  CHK_ERROR(err);
  err = clWaitForEvents(1, &event);
  CHK_ERROR(err);
  err = clFinish(cmd_queue);
  CHK_ERROR(err);

  cl_ulong time_start;
  cl_ulong time_end;
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START,
      sizeof(time_start), &time_start, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,
      sizeof(time_end), &time_end, NULL);
  double nano_seconds = time_end - time_start;

#ifdef ARRAY_SIZE
	printf("Done in %f seconds.\n", nano_seconds / 1e9);
#else
	printf("%f\n", nano_seconds / 1e9);
#endif

  // Finally, release all that we have allocated.
  err = clReleaseCommandQueue(cmd_queue);
  CHK_ERROR(err);
  err = clReleaseContext(context);
  CHK_ERROR(err);
  free(platforms);
  free(device_list);

  /* Check result */
  if (!assert_equal(Y_reference, Y, ARRAY_SIZE)) {
    printf("Error: arrays are different\n");
  } else {
#ifdef ARRAY_SIZE
    printf("Success: arrays are equal\n");
#endif
  }

  
  return 0;
}



// The source for this particular version is from: https://stackoverflow.com/questions/24326432/convenient-way-to-show-opencl-error-codes
const char* clGetErrorString(int errorCode) {
  switch (errorCode) {
  case 0: return "CL_SUCCESS";
  case -1: return "CL_DEVICE_NOT_FOUND";
  case -2: return "CL_DEVICE_NOT_AVAILABLE";
  case -3: return "CL_COMPILER_NOT_AVAILABLE";
  case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
  case -5: return "CL_OUT_OF_RESOURCES";
  case -6: return "CL_OUT_OF_HOST_MEMORY";
  case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
  case -8: return "CL_MEM_COPY_OVERLAP";
  case -9: return "CL_IMAGE_FORMAT_MISMATCH";
  case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
  case -12: return "CL_MAP_FAILURE";
  case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
  case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
  case -15: return "CL_COMPILE_PROGRAM_FAILURE";
  case -16: return "CL_LINKER_NOT_AVAILABLE";
  case -17: return "CL_LINK_PROGRAM_FAILURE";
  case -18: return "CL_DEVICE_PARTITION_FAILED";
  case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
  case -30: return "CL_INVALID_VALUE";
  case -31: return "CL_INVALID_DEVICE_TYPE";
  case -32: return "CL_INVALID_PLATFORM";
  case -33: return "CL_INVALID_DEVICE";
  case -34: return "CL_INVALID_CONTEXT";
  case -35: return "CL_INVALID_QUEUE_PROPERTIES";
  case -36: return "CL_INVALID_COMMAND_QUEUE";
  case -37: return "CL_INVALID_HOST_PTR";
  case -38: return "CL_INVALID_MEM_OBJECT";
  case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
  case -40: return "CL_INVALID_IMAGE_SIZE";
  case -41: return "CL_INVALID_SAMPLER";
  case -42: return "CL_INVALID_BINARY";
  case -43: return "CL_INVALID_BUILD_OPTIONS";
  case -44: return "CL_INVALID_PROGRAM";
  case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
  case -46: return "CL_INVALID_KERNEL_NAME";
  case -47: return "CL_INVALID_KERNEL_DEFINITION";
  case -48: return "CL_INVALID_KERNEL";
  case -49: return "CL_INVALID_ARG_INDEX";
  case -50: return "CL_INVALID_ARG_VALUE";
  case -51: return "CL_INVALID_ARG_SIZE";
  case -52: return "CL_INVALID_KERNEL_ARGS";
  case -53: return "CL_INVALID_WORK_DIMENSION";
  case -54: return "CL_INVALID_WORK_GROUP_SIZE";
  case -55: return "CL_INVALID_WORK_ITEM_SIZE";
  case -56: return "CL_INVALID_GLOBAL_OFFSET";
  case -57: return "CL_INVALID_EVENT_WAIT_LIST";
  case -58: return "CL_INVALID_EVENT";
  case -59: return "CL_INVALID_OPERATION";
  case -60: return "CL_INVALID_GL_OBJECT";
  case -61: return "CL_INVALID_BUFFER_SIZE";
  case -62: return "CL_INVALID_MIP_LEVEL";
  case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
  case -64: return "CL_INVALID_PROPERTY";
  case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
  case -66: return "CL_INVALID_COMPILER_OPTIONS";
  case -67: return "CL_INVALID_LINKER_OPTIONS";
  case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";
  case -69: return "CL_INVALID_PIPE_SIZE";
  case -70: return "CL_INVALID_DEVICE_QUEUE";
  case -71: return "CL_INVALID_SPEC_ID";
  case -72: return "CL_MAX_SIZE_RESTRICTION_EXCEEDED";
  case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
  case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
  case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
  case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
  case -1006: return "CL_INVALID_D3D11_DEVICE_KHR";
  case -1007: return "CL_INVALID_D3D11_RESOURCE_KHR";
  case -1008: return "CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR";
  case -1009: return "CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR";
  case -1010: return "CL_INVALID_DX9_MEDIA_ADAPTER_KHR";
  case -1011: return "CL_INVALID_DX9_MEDIA_SURFACE_KHR";
  case -1012: return "CL_DX9_MEDIA_SURFACE_ALREADY_ACQUIRED_KHR";
  case -1013: return "CL_DX9_MEDIA_SURFACE_NOT_ACQUIRED_KHR";
  case -1093: return "CL_INVALID_EGL_OBJECT_KHR";
  case -1092: return "CL_EGL_RESOURCE_NOT_ACQUIRED_KHR";
  case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
  case -1057: return "CL_DEVICE_PARTITION_FAILED_EXT";
  case -1058: return "CL_INVALID_PARTITION_COUNT_EXT";
  case -1059: return "CL_INVALID_PARTITION_NAME_EXT";
  case -1094: return "CL_INVALID_ACCELERATOR_INTEL";
  case -1095: return "CL_INVALID_ACCELERATOR_TYPE_INTEL";
  case -1096: return "CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL";
  case -1097: return "CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL";
  case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
  case -1098: return "CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL";
  case -1099: return "CL_INVALID_VA_API_MEDIA_SURFACE_INTEL";
  case -1100: return "CL_VA_API_MEDIA_SURFACE_ALREADY_ACQUIRED_INTEL";
  case -1101: return "CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL";
  default: return "CL_UNKNOWN_ERROR";
  }
}
