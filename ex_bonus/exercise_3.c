// Template file for the OpenCL Assignment 4

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
#include <sys/time.h>

// Number of dimensions used to specify work-items in a work-group
#define WORK_DIM 			1
#define KERNEL_FILEPATH 	"particle_simulator.cl"
#define KERNEL_NAME			"particle_simulator"

#define ABS(a) ((a) < 0 ? -(a) : (a))

// This is a macro for checking the error variable.
#define CHK_ERROR(err) if (err != CL_SUCCESS) fprintf(stderr,"Error: %s\n",clGetErrorString(err));

// Number of particles to simulate
int NUM_PARTICLES = 16384;

// Number of times to simulate the movement (over 1 second) of a particle
int NUM_ITERATIONS = 1;

// Number of work-items per work-group
int BLOCK_SIZE = 128;

// Gravity field for updating particle velocities
cl_float3 field = (cl_float3) {0.f, 0.f, 9.8f};

// Structure for particles
typedef struct {
	cl_float3 position;
	cl_float3 velocity;
} Particle;

// A errorCode to string converter (forward declaration)
const char* clGetErrorString(int);

/**
 * Change the position of the given particle based on its velocity using the
 * formula `new_position.coord = old_position.coord + velocity.coord` where
 * coord is x, y and z.
 *
 * @param particle	Particle for which a position update will be performed
 */
void updatePosition(Particle *particle) {
	particle->position.x = particle->position.x + particle->velocity.x;
	particle->position.y = particle->position.y + particle->velocity.y;
	particle->position.z = particle->position.z + particle->velocity.z;
}

/**
 * Update the velocity of the given particle according to a field that specifies
 * the rate of change for each dimension of the particle's velocity
 *
 * @param particle	Particle for which a velocity update will be performed
 * @param field		Rate of change for each dimension (x, y, z) of a velocity
 */
void updateVelocity(Particle *particle, cl_float3 field) {
	particle->velocity.x = particle->velocity.x + field.x;
	particle->velocity.y = particle->velocity.y + field.y;
	particle->velocity.z = particle->velocity.z + field.z;
}

/**
 * Host implementation for the simulation of moving particles
 *
 * @param particles			List of particles for which to simulate movement
 * @param num_particles		# of particles to simulate
 * @param num_iterations	# of timesteps for which to simulate each particle
 */
void simulateParticlesHost(Particle *particles,
    int num_particles, int num_iterations) {
  
	for (Particle *particle = particles;
		particle < particles + num_particles;
		particle++) {
    	for (int i = 0; i < num_iterations; ++i) {
			// Update velocity first
			updateVelocity(particle, field);

			// Update position
			updatePosition(particle);
		}
	}
}

/**
 * Fill the given array with n random floats.
 *
 * @param array	Array to populate with floats.
 * @param n		Number of floats to populate the array with.
 */
void populateParticleArray(Particle *particles, int n) {
	Particle particle;

	for (int index = 0; index < n; index++) {
		// Generate random particles
		particle.position.x = 10.0 * ((float) rand() / (float) RAND_MAX);
		particle.position.y = 10.0 * ((float) rand() / (float) RAND_MAX);
		particle.position.z = 10.0 * ((float) rand() / (float) RAND_MAX);
		particle.velocity.x = 1.0 * ((float) rand() / (float) RAND_MAX);
		particle.velocity.y = 1.0 * ((float) rand() / (float) RAND_MAX);
		particle.velocity.z = 1.0 * ((float) rand() / (float) RAND_MAX);
		particles[index] = particle;
	}
}

/**
 * Compare the simulation results of the device and host implementations.
 *
 * @param deviceOut	Outcome of simulation from the device implementation
 * @param hostOut	Outcome of simulation from the host implementation
 * @param n			The size of deviceOut and hostOut
 */
void compareSimulationResults(Particle *deviceOut, Particle *hostOut, int n) {
	int resultsAreEqual = 1;
	printf("Comparing the output for each implementation... ");

	for (int index = 0; index < n; index++) {
		float cumDiff = 0;
		cumDiff += ABS(deviceOut[index].position.x - hostOut[index].position.x);
		cumDiff += ABS(deviceOut[index].position.y - hostOut[index].position.y);
		cumDiff += ABS(deviceOut[index].position.z - hostOut[index].position.z);
		cumDiff += ABS(deviceOut[index].velocity.x - hostOut[index].velocity.x);
		cumDiff += ABS(deviceOut[index].velocity.y - hostOut[index].velocity.y);
		cumDiff += ABS(deviceOut[index].velocity.z - hostOut[index].velocity.z);

		// Difference is larger than rounding-error tolerance of .001, means
		// the outcomes are too different
		if (cumDiff > .001 || cumDiff < -.001) {
			resultsAreEqual = 0;
			break;
		}
	}

	// The outcomes of simulation for the device and host implementations are equal
	if (resultsAreEqual) {
		printf("Correct!\n");
	} else {
		printf("INCORRECT!!!\n");
	}
}

/**
 * Given the filepath for a kernel, open, read and return the source code for
 * the kernel.
 *
 * @param kernelFilepath	Filepath of the kernel to return source code for
 *
 * @return the source code (as a character array) of the kernel
 */
char* getKernelSourceCode(char *kernelFilepath) {
	long fileLength;
	char *sourceCode = 0;
	FILE *kernelFile = fopen(kernelFilepath, "rb");

	// File does not exist or is not accessible
	if (kernelFile == NULL) {
		fprintf(stderr, "Could not open %s\n", kernelFilepath);
		return NULL;
	}

	// Read to the end of the file in order to retrieve the length of the file
	// using ftell
	fseek(kernelFile, 0, SEEK_END);
	fileLength = ftell(kernelFile);
	
	// Move back to the beginning of the file for reading the data
	fseek(kernelFile, 0, SEEK_SET);

	// Create enough memory to hold the entire file
	sourceCode = (char *) malloc(fileLength);

	// Could not allocate space for sourceCode
	if (sourceCode == NULL) {
		fprintf(stderr, "Malloc failed\n");
		return NULL;
	}

	// Read the file into the sourceCode buffer
	fread (sourceCode, 1, fileLength, kernelFile);
	fclose (kernelFile);
	return sourceCode;

}

/**
 * Return a timestamp with double precision
 */
double cpuSecond() {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

/**
 * Execute and time the host implementation of particle simulation.
 */
void runHostSimulation(Particle *hostParticles) {
	// Execute and time the host particle simulator
	double startTime = cpuSecond();
	simulateParticlesHost(hostParticles, NUM_PARTICLES, NUM_ITERATIONS);
	printf("CPU implementation took %f seconds\n", cpuSecond() - startTime);
}

void runDeviceSimulation(Particle *hostParticles, Particle *saveDevParticles,
		int particleArraySize) {
	/* BEGIN OPENCL OVERHEAD */
	cl_platform_id * platforms; cl_uint n_platform;

	// Find OpenCL Platforms
	cl_int err = clGetPlatformIDs(0, NULL, &n_platform); CHK_ERROR(err);
	platforms = (cl_platform_id *) malloc(sizeof(cl_platform_id)*n_platform);
	err = clGetPlatformIDs(n_platform, platforms, NULL); CHK_ERROR(err);

	// Find and sort devices
	cl_device_id *device_list; cl_uint n_devices;
	err = clGetDeviceIDs( platforms[0], CL_DEVICE_TYPE_GPU, 0,NULL, &n_devices);
	CHK_ERROR(err);
	device_list = (cl_device_id *) malloc(sizeof(cl_device_id)*n_devices);
	err = clGetDeviceIDs(
		platforms[0], CL_DEVICE_TYPE_GPU, n_devices, device_list, NULL);
	CHK_ERROR(err);

	// Create and initialize an OpenCL context
	cl_context context = clCreateContext(
		NULL, n_devices, device_list, NULL, NULL, &err); CHK_ERROR(err);

	// Create a command queue
	cl_command_queue cmd_queue = clCreateCommandQueue(
			context, device_list[0], 0, &err); CHK_ERROR(err); 

	// Fetch the source code of the kernel found at KERNEL_FILEPATH
	char* sourceCode = getKernelSourceCode(KERNEL_FILEPATH);

	// Create an OpenCL program
	cl_program program = clCreateProgramWithSource(
		context, 1, (const char**) &sourceCode, NULL, &err); CHK_ERROR(err);

	err = clBuildProgram(program, 1, device_list, NULL, NULL, NULL);
	CHK_ERROR(err);

	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];
		clGetProgramBuildInfo(program, device_list[0],
			CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("What da heck\n");
		fprintf(stderr, "Build error: %s\n", buffer);
		printf("seriously\n");
		printf("%s\n", buffer);
		exit(0);
	}

	// Create a kernel object referencing our KERNEL_NAME kernel
	cl_kernel kernel = clCreateKernel(program, KERNEL_NAME, &err);
	CHK_ERROR(err);

	/* END OPENCL OVERHEAD */

	/* BEGIN EXECUTION OF THE DEVICE IMPLEMENTATION */
	// Allocate device memory
	cl_mem deviceParticles = clCreateBuffer(context, CL_MEM_READ_WRITE,
			particleArraySize, NULL, &err); CHK_ERROR(err);
	
	// Write to device memory
	err = clEnqueueWriteBuffer(cmd_queue, deviceParticles, CL_TRUE, 0,
			particleArraySize, hostParticles, 0, NULL, NULL);
	CHK_ERROR(err);

	// Set arguments for the kernel
	CHK_ERROR(
		clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &deviceParticles));
	CHK_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_float3), (void *) &field));
	CHK_ERROR(clSetKernelArg(kernel, 2, sizeof(int), (void *) &NUM_PARTICLES));
	CHK_ERROR(clSetKernelArg(kernel, 3, sizeof(int), (void *) &NUM_ITERATIONS));

	/* Parameters about number of work groups and number of work items per work
	* group to execute the kernel with */
	size_t n_workitem[WORK_DIM] = {NUM_PARTICLES};
	size_t workgroup_size[WORK_DIM] = {BLOCK_SIZE};

	// Launch the kernel
	err = clEnqueueNDRangeKernel(cmd_queue, kernel, WORK_DIM, NULL, n_workitem,
			workgroup_size, 0, NULL, NULL); CHK_ERROR(err);

	// Get result back from the host
	err = clEnqueueReadBuffer(cmd_queue, deviceParticles, CL_TRUE, 0,
			particleArraySize, saveDevParticles, 0, NULL, NULL);
	CHK_ERROR(err);

	/* Wait and make sure everything finished */
	err = clFlush(cmd_queue); CHK_ERROR(err);
	err = clFinish(cmd_queue); CHK_ERROR(err);
	/* END EXECUTION OF THE DEVICE IMPLEMENTATION */

	// Finally, release all that we have allocated.
	err = clReleaseCommandQueue(cmd_queue); CHK_ERROR(err);
	err = clReleaseContext(context); CHK_ERROR(err);
	err = clReleaseProgram(program); CHK_ERROR(err);
	err = clReleaseKernel(kernel); CHK_ERROR(err);
	free(platforms);
	free(device_list);
}

void parseArguments(int argc, char **argv) {
	NUM_PARTICLES 	= argc > 1 ? atoi(argv[1]) : NUM_PARTICLES;
	NUM_ITERATIONS 	= argc > 2 ? atoi(argv[2]) : NUM_ITERATIONS;
	BLOCK_SIZE 	= argc > 3 ? atoi(argv[3]) : BLOCK_SIZE;

	printf("Running simulation with:\n\t"
			"%d particles\n\t"
			"%d iterations\n\t"
			"%d block size\n",
			NUM_PARTICLES, NUM_ITERATIONS, BLOCK_SIZE
	);
}

int main(int argc, char **argv) {
	// Get values for number of particles, iteration, and block size, or quit
	parseArguments(argc, argv);

	int particleArraySize = sizeof(Particle) * NUM_PARTICLES;

	// Allocate space to hold host particles and then populate the space
	Particle *hostParticles = (Particle *) malloc(particleArraySize);

	populateParticleArray(hostParticles, NUM_PARTICLES);

	// Allocate space to hold the result of simulating device particles
	Particle *devResult = (Particle *) malloc(particleArraySize);
	
	// Measure how long it takes to run the device implementation
	runDeviceSimulation(hostParticles, devResult, particleArraySize);
	
	// Measure how long it takes to run the host implementation
	runHostSimulation(hostParticles);

	// Compare the results of the device and host simulations
	compareSimulationResults(devResult, hostParticles, particleArraySize);

	// Free all memory
	free(hostParticles);
	free(devResult);

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
