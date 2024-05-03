# Use an official C++ environment as the base image
FROM gcc:latest

# Install required packages
RUN apt-get update && apt-get install -y \
    cmake \
    libeigen3-dev \
    libomp-dev \
    git \
 && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /usr/src/clipperplus

# Copy your project into the container
COPY . .

# Configure the project and enable tests
#RUN cmake -S . -B build -DBUILD_TESTS=ON

# Build the project and the tests
#RUN cmake --build build

# Run the tests
#RUN cd build && ctest || true

#CMD ["./build/test/cpp_tests/tests"]
#CMD ["./build/clipperplus_test"]
CMD ["bash"]
