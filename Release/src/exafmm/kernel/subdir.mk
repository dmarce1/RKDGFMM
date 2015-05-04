################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../src/exafmm/kernel/CPUSphericalLaplace.cxx 

OBJS += \
./src/exafmm/kernel/CPUSphericalLaplace.o 

CXX_DEPS += \
./src/exafmm/kernel/CPUSphericalLaplace.d 


# Each subdirectory must supply rules for building sources it contributes
src/exafmm/kernel/%.o: ../src/exafmm/kernel/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -march=corei7-avx -std=c++0x -DNDEBUG -DSpherical -DEXPANSION=5 -I"/home/dmarce1/workspace/RKDGFMM/src/exafmm/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


