# syntax=docker/dockerfile:1

# Build imgage for Trilinos Sacado
# This will be resued by https://gitlab.com/bootstrapcollaboration/scalar_blocks

# latest alpine release
FROM alpine:3.18 AS build

RUN apk add \
    cmake \
    g++ \
    gfortran \
    git \
    make \
    python3
WORKDIR /usr/local/src/Trilinos
# Build Trilinos Sacado from current sources
COPY . .
RUN mkdir -p build && \
    cd build && \
    cmake .. -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Kokkos=OFF -DTrilinos_ENABLE_Teuchos=OFF \
          -DTrilinos_ENABLE_TESTS=ON \
          -DCMAKE_INSTALL_PREFIX=/usr/local/ && \
    make && \
    make install

# Take only Trilinos headers and binaries + load necessary dynamic libraries
FROM alpine:3.18 as install
RUN apk add libstdc++
COPY --from=build /usr/local/bin /usr/local/bin
COPY --from=build /usr/local/include /usr/local/include
COPY --from=build /usr/local/lib /usr/local/lib

# Test target, usage:
# docker build . --tag trilinos-sacado-test --target test
# docker run trilinos-sacado-test ctest --output-on-failure --no-tests=error
FROM install as test
RUN apk add cmake # for ctest
COPY --from=build /usr/local/src/Trilinos/build /usr/local/src/Trilinos/build
WORKDIR /usr/local/src/Trilinos/build

FROM install