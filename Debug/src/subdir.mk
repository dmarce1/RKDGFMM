################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/conserved.cpp \
../src/exafmm.cpp \
../src/fourier_legendre.cpp \
../src/grid.cpp \
../src/initial.cpp \
../src/legendre.cpp \
../src/main.cpp \
../src/primitive.cpp 

OBJS += \
./src/conserved.o \
./src/exafmm.o \
./src/fourier_legendre.o \
./src/grid.o \
./src/initial.o \
./src/legendre.o \
./src/main.o \
./src/primitive.o 

CPP_DEPS += \
./src/conserved.d \
./src/exafmm.d \
./src/fourier_legendre.d \
./src/grid.d \
./src/initial.d \
./src/legendre.d \
./src/main.d \
./src/primitive.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -D_GLIBCXX_DEBUG -I/home/dmarce1/include -I/usr/include/mpi -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


