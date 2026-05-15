# Chromosphere Model Library

This a stream-aligned chromosphere model in C++. The model uses Armadillo, a high-performance linear algebra library, and provides various utilities and functions for simulating chromospheric phenomena.


## Naming Standard

#### Vec (arma::col<float>)

1. 1-D vectors contain `ns` elements,  and are named as `A_i(ph/mh)`, where `_i` denotes it has 1 dimension.
2. 2-D matrices are flattened into 1-D Vec vectors with the size of `ns*num_of_eq`. They are named as `A_ii(ph/mh)`.
3. For vectors without postfix `_i` or `_ii`, the dimension is uncertain, but can be determined by the context in most cases without confusion.

## Table of Contents

- [Chromosphere Namespace](#chromosphere-namespace)
  - [Constants](#constants)
  - [Variables](#variables)
  - [Functions](#functions)
  - [RHS Solver](#rhs-solver)
- [Other Utilities](#other-utilities)
- [Fortran Interfaces](#fortran-interfaces)
- [ModChromosphereCPP](#ModChromosphereCPP)

## Chromosphere Namespace

The `chromosphere` namespace encapsulates all the constants, variables, and functions used in the chromosphere model. It is designed to handle the specific needs of modeling chromospheric processes using finite difference methods.

### Constants

These constants are predefined within the namespace and are used throughout the simulation.

- **Number of Equations**: `const uword num_of_eq = 6;`
- **Finite Difference Methods**:
  - `FORWARD = 0`
  - `BACKWARD = 1`
  - `CENTRAL = 2`
- **Permutation Functions**:
  - `SLICE = 1`
  - `CUBE = num_of_eq`
- **Small Enough Threshold**: `const float SMALL_ENOUGH = 1.0e-6;`

### Variables

These variables are extern and are defined elsewhere in the implementation.

- **Grid Parameters**:
  - `uword ns;` - Number of grid points in the x-direction
  - `uword num_of_elem;` - Total number of elements (ns * num_of_eq)

- **Physical Constants**:
  - `float gammamono;` - Ratio of specific heats
  - `float alpha_p;` - Divergence error propagation parameter
  - `float pi;` - The value of π (pi)
  - `float m_i;` - Mass of proton in kilograms
  - `float mu_0;` - Permeability of free space
  - `float k_b;` - Boltzmann constant in SI units
  - `float e_;` - Elementary charge in Coulombs

- **Domain Normalization Constants**:
  - `float dom_l_0;` - Domain length normalization constant in meters
  - `float dom_t_0;` - Time normalization constant in seconds
  - `float width;` - Width of the domain in normalized units
  - `float height;` - Height of the domain in normalized units
  - `Vec ds_i;` - Length of i-th cell
  - `float CFL;` - Courant-Friedrichs-Lewy condition

- **Normalized Physical Quantities**:
  - `float n_0;` - Density normalization constant in m^{-3}
  - `float B_0;` - Magnetic field normalization constant in Teslas
  - `float t_0;` - Time normalization constant in seconds
  - `float l_0;` - Length normalization constant in meters
  - `float m_n;` - Normalized mass of neutrals
  - `float m_e;` - Normalized mass of electron
  - `float V_0;` - Alfven velocity normalization constant in m/s
  - `float p_0;` - Pressure normalization constant in Pascals
  - `float Temp_0;` - Temperature normalization constant in Kelvin
  - `float g;` - Gravitational acceleration in normalized units
  - `float q_0;` - Heat transfer normalization constant
  - `float e_0;` - Charge normalization constant in Coulombs
  - `float q;` - Normalized elementary charge constant
  - `float kappa_0;` - Heat transfer coefficient normalization constant

- **Stream-Aligned Quantities**:
  - `Vec cellX_i, cellY_i, cellZ_i;` - Cell coordinates XYZ.
  - `Vec B_imh, B_iph, B_i;` - Magnetic field values at various cell faces.
  - `Vec dinvB_ds_i;` - Derivative of 1/B at i-th cell.

- **Broadcasted Variables**:
  - `Vec B_iimh, B_iiph, B_ii, dinvB_ds_ii, dt_ii, ds_ii;` - Precomputed variables to avoid redundant calculations.

### Functions

The functions within the `chromosphere` namespace handle the initialization, computation, and advancement of the chromosphere model.

- **Initialization**:
  - `void init_params(uword input_nx);`
  - `Vec shock_tube_ic();` - Initialize conditions for shock tube problem.

- **Advancement Methods**:
  - `Vec cal_dt_i(const Vec& xn);` - Calculate local time steps.
  - `Vec advance_Euler_ii(const Vec& xn, const Vec& dt);` - Euler method for time advancement.
  - `Vec advance_RK4(const Vec& xn, const float& dt);` - RK4 method for time advancement.

- **Collision Rates**:
  - `float nu_nn(const float& nn, const float& T);`
  - `Vec nu_nn(const Vec& nn, const Vec& T);`
  - `float nu_in(const float& nn, const float& T);`
  - `Vec nu_in(const Vec& nn, const Vec& T);`
  - `Vec nu_en(const Vec& nn, const Vec& T);` - (To be completed)
  - `Vec nu_ei(const Vec& ni, const Vec& Te);`

- **Permutation Functions**:
  - `Vec ip1(const Vec& xn, const uword nk = CUBE);` - Permute the vector `xn` one position to the right.
  - `Vec im1(const Vec& xn, const uword nk = CUBE);` - Permute the vector `xn` one position to the left.

- **Derivatives and Norms**:
  - `Vec deriv_s(const Vec& xn, const uword type = FORWARD);` - Compute the derivative of `xn` using finite difference methods (forward, backward, or central).
  - `Vec vec_abs(const Vec& Ax, const Vec& Ay, const Vec& Az);` - Calculate the magnitude of a vector given its components.
  - `Vec vec_abs2(const Vec& Ax, const Vec& Ay, const Vec& Az);` - Calculate the squared magnitude of a vector.

- **Conversions**:
  - `Vec cons2prim(const Vec& cons);` - Convert conserved variables to primitive variables.
  - `Vec prim2cons(const Vec& prim);` - Convert primitive variables to conserved variables.

### RHS Solver

The RHS (Right-Hand Side) Solver functions are used to compute the explicit and implicit terms in the governing equations of the chromosphere model.

- `Vec flux_lim(const Vec& r);` - Apply a flux limiter to the vector `r` to prevent numerical oscillations.
- `Vec rhs_explicit_ii(const Vec& xn_ii, const float& dt_i);` - Calculate the explicit terms in the right-hand side of the equation for a given time step.
- `Vec rhs_implicit_ii(const Vec& xn_ii);` - Calculate the implicit terms in the right-hand side of the equation.
- `Vec find_spectral_radius_ii(const Vec& xn_ii);` - Calculate the fast speed eigenvalue (spectral radius) in the x-direction.
- `Vec cal_F_ii(const Vec& xn_ii);` - Compute the flux in the x-direction for the given state vector `xn_ii`.
- `Vec get_max_v_i(const Vec& xn_ii);` - Get the maximum velocity for Courant condition calculation from the state vector `xn_ii`.

## Other Utilities

The header file includes additional utility functions for converting between arrays and vectors:

- `void array2vec(const float* from_array, Vec& to_vec, const uword& size);` - Convert a C-style array to an Armadillo vector.
- `void vec2array(const Vec& from_vec, float* to_array, const uword& size);` - Convert an Armadillo vector to a C-style array.

## Fortran Interfaces

To facilitate interaction with Fortran code, several interfaces are provided. Below is a description of each function, including the intents of their arguments:

- `extern "C" void chromo_state_();`
  - **Description**: Gathers the state for the Space Weather Modeling Framework (SWMF).
  - **Arguments**: None

- `extern "C" int Chromo_sub2ind(const int* i, const int* k);`
  - **Description**: Converts a subscript pair `(i, k)` to a flattened vector index for accessing a 2D array stored in a 1D array.
  - **Arguments**:
    - `i` (intent: in) - Row index.
    - `k` (intent: in) - Column index.

- `extern "C" void Chromo_init_CPP(float* xn, const int* ns, struct Chromo_params* params);`
  - **Description**: Initializes the parameters needed for the simulation, optionally using the `xn` and `params` structures.
  - **Arguments**:
    - `xn` (intent: inout) - Array containing the initial state variables.
    - `ns` (intent: in) - Number of grid points in the x-direction.
    - `params` (intent: inout) - Structure containing the parameters for the simulation.

- `extern "C" void Chromo_cal_dt(const float* xn, float* dt);`
  - **Description**: Calculates local time stepping `dt` for each grid cell.
  - **Arguments**:
    - `xn` (intent: in) - Array containing the current state variables.
    - `dt` (intent: out) - Array to store the calculated time steps for each cell.

- `extern "C" void Chromo_advance_Euler(float* xn, const float* dt);`
  - **Description**: Advances the simulation state using the Euler method for time integration.
  - **Arguments**:
    - `xn` (intent: inout) - Array containing the current state variables, which will be updated.
    - `dt` (intent: in) - Array containing the time step for each cell.

- `extern "C" void Chromo_prim2cons(float* p_xn);`
  - **Description**: Converts primitive variables to conserved variables.
  - **Arguments**:
    - `p_xn` (intent: inout) - Array of primitive variables to be converted to conserved variables.

- `extern "C" void Chromo_cons2prim(float* c_xn);`
  - **Description**: Converts conserved variables to primitive variables.
  - **Arguments**:
    - `c_xn` (intent: inout) - Array of conserved variables to be converted to primitive variables.

## ModChromosphereCPP

The `chromo_params_F` type in the `ModChromosphereCPP` module is a Fortran data structure designed to represent the physical parameters, cell parameters, and state variables of the chromosphere model in a way that is compatible with Fortran-style arrays and operations. This structure is used in the Fortran portion of the code, while the equivalent `chromo_params_C` type is used in the C++ code. The module includes routines to convert between these two representations.

### Structure of `chromo_params_F`

- **Physical Parameters**: These fields represent various physical constants and normalization factors used in the model.
  - `gammamono`, `alpha_p`, `pi`: Adiabatic index, divergence error propagation parameter, and the value of π.
  - `m_i`, `mu_0`, `k_b`, `e_`: Mass of the proton, permeability of free space, Boltzmann constant, and elementary charge.
  - `dom_l_0`, `dom_t_0`, `width`, `height`: Domain normalization constants for length, time, and dimensions of the domain.
  - `n_0`, `B_0`, `t_0`, `l_0`, `m_n`, `m_e`: Normalization constants for density, magnetic field, time, length, neutral mass, and electron mass.
  - `V_0`, `p_0`, `Temp_0`, `g`, `q_0`, `e_0`, `q`, `kappa_0`: Normalization constants for Alfven velocity, pressure, temperature, gravity, heat transfer, charge, and heat transfer coefficient.

- **Cell Parameters**: These fields store indices and parameters related to the computational grid and the conservation equations.
  - `CNI`, `CNN`, `CNV`, `CNU`, `CEI`, `CEN`: Indices for conserved variables, including ion and neutral densities, velocities, and energies.
  - `PNI`, `PNN`, `PV`, `PU`, `PPI`, `PPN`: Indices for primitive variables, including ion and neutral densities, velocities, and pressures.
  - `ns`, `num_of_eq`, `num_of_elem`: Number of grid points in the x-direction, number of equations, and total number of elements.
  - `CFL`: Courant-Friedrichs-Lewy condition for time-stepping stability.

- **Dynamic Arrays**: These pointers reference arrays that store the grid spacing, cell positions, magnetic field values, and state variables.
  - `LengthSi_G`: Pointer to an array representing the grid spacing in the x-direction.
  - `Coord_DF`: Pointer to a 2D array representing cell positions in the domain (X, Y, Z).
  - `BSi_F`: Pointer to an array storing magnetic field values.
  - `State_VC`: Pointer to a 2D array representing the state variables (conserved or primitive).
  - `OuterBoundary`: Pointer to a 2D array representing boundary conditions at the outer edges of the chromosphere.

- **Memory Management**:
  - `is_allocated`: A logical flag indicating whether the dynamic arrays have been allocated.
  - The module provides subroutines `Chromo_allocate_params_F` and `Chromo_deallocate_params_F` to manage the memory allocation and deallocation of these arrays.

### Key Routines

- **Chromo_ConvertParam_C2F**: This subroutine converts data from the `chromo_params_C` structure (used in C++) to the `chromo_params_F` structure (used in Fortran). It handles the conversion of C-style pointers to Fortran arrays using the `c_f_pointer` intrinsic, which allows seamless integration between C++ and Fortran code.

- **Chromo_ConvertParam_F2C**: Although the detailed implementation is not provided, this subroutine is intended to convert data from the `chromo_params_F` structure back to the `chromo_params_C` structure.

- **Chromo_allocate_params_F**: This subroutine allocates the necessary memory for the dynamic arrays within `chromo_params_F`, ensuring that all required arrays are properly initialized before use.

- **Chromo_deallocate_params_F**: This subroutine deallocates the memory associated with the dynamic arrays within `chromo_params_F`, freeing up resources when they are no longer needed.

- **Chromo_init**: This subroutine initializes the chromosphere model by converting parameters between the Fortran and C++ structures, calling the C++ initialization routine (`Chromo_init_CPP`), and then converting the parameters back for continued use in Fortran.

## Data Transfer

1. **Fortran-side Storage**: `chromo_params_F` holds the data on the Fortran side.
2. **Data Transfer**: The function `Chromo_ConvertParam_F2C` is responsible for transferring data from `chromo_params_F` to `chromo_params_C`. Note that `xn_I` has a dimension of (ns, n_eq) while `State_CV` is (n_eq, ns).
3. **C++ Binding**: The `chromo_params_C` data is bound to the `struct Chromo_params` in `chromo_interface.cpp`, allowing it to be passed as an argument within the C++ code.
4. **Solver Integration**: Within `chromo_interface.cpp`, the contents of `struct Chromo_params` are copied into variables within the `namespace chromosphere`, which are then directly utilized by the solver.

Notes:

When transferring data, we use two subroutines: `c_f_pointer` and `c_loc`. They associate C pointers and fortran arrays without allocating new memory. In short, if you have a `c_ptr` variable `A`  and you want to access its elements in Fortran, you use `c_f_pointer(A,B,size)` to assign the address to a Fortran array `B`; if you have a Fortran array `B` and you want to transform it into a `c_ptr` type so it can be accessed by C code, you use `A=c_loc(B)`. As a result, `A` and `B` point to the same address in the memory.

It's worth noting that the interoperability is not guaranteed when you try to pass arrays. So the following code may not work.

```fortran
subroutine Chromo_cal_dt_CPP(pC, dt_I) bind(C, name="Chromo_cal_dt")
  use iso_c_binding,     ONLY: c_float
  import chromo_params_C
  type(chromo_params_C), intent(in) :: pC
  real(c_float), dimension(:), intent(out) :: dt_I
end subroutine Chromo_cal_dt_CPP
```

and

```c++
extern "C" void Chromo_cal_dt(const struct Chromo_params* params, float* dt);
```

The correct way is

```fortran
subroutine Chromo_cal_dt_CPP(pC, dt_I) bind(C, name="Chromo_cal_dt")
  use iso_c_binding,     ONLY: c_float, c_ptr
  import chromo_params_C
  type(chromo_params_C), intent(in) :: pC
  type(c_ptr), value, intent(out) :: dt_I
end subroutine Chromo_cal_dt_CPP
```

Here the `value` attribute is used bacause otherwise Fortran will pass the address of `type(c_ptr)` to `Chromo_cal_dt_CPP`, which, to C, means `float **`.



If you do not believe it, try this example. The goal is to exchange an array between Fortran and C++ and modify it in both environments.

```fortran
! FORTRAN_PART.F90
MODULE InteropModule
    USE ISO_C_BINDING
    IMPLICIT NONE

    TYPE, BIND(C) :: ArrayStruct
        TYPE(c_ptr) :: data_ptr   ! Pointer to an integer array
        INTEGER(C_INT) :: length  ! Size of the array
    END TYPE ArrayStruct

    INTERFACE
        SUBROUTINE ProcessArrayInC(arr_struct) BIND(C, name='C_subroutine')
            USE ISO_C_BINDING
            IMPORT ArrayStruct
            TYPE(ArrayStruct), INTENT(INOUT) :: arr_struct
        END SUBROUTINE ProcessArrayInC
    END INTERFACE

END MODULE InteropModule

PROGRAM MainProgram
    USE InteropModule
    USE ISO_C_BINDING, ONLY: c_loc, c_f_pointer, c_int
    IMPLICIT NONE

    TYPE(ArrayStruct) :: array_struct
    INTEGER(c_int), POINTER :: fortran_array(:)
    INTEGER(c_int) :: i

    ! Allocate and initialize the Fortran array
    ALLOCATE(fortran_array(5))
    FORALL(i = 1:5) fortran_array(i) = 5 - i

    ! Set the array size and point the C pointer to the Fortran array
    array_struct%length = 5
    array_struct%data_ptr = c_loc(fortran_array)

    ! Pass the array to the C++ subroutine for processing
    CALL ProcessArrayInC(array_struct)

    ! Convert the C pointer back to a Fortran pointer and print the modified array
    CALL c_f_pointer(array_struct%data_ptr, fortran_array, [5])
    PRINT *, 'Modified array after C++ subroutine call:', fortran_array

END PROGRAM MainProgram
```

```C++
// CPP_PART.CPP
#include <iostream>
using namespace std;

struct ArrayStruct {
    int* data_ptr;  // Pointer to an integer array
    int length;     // Size of the array
};

extern "C" {
    void C_subroutine(ArrayStruct* arr_struct);
}

extern "C" void C_subroutine(ArrayStruct* arr_struct) {
    // Modify the array elements in C++
    for (int i = 0; i < arr_struct->length; ++i) {
        cout << "Element " << i << ": " << arr_struct->data_ptr[i] << endl;
        arr_struct->data_ptr[i] += 10;
    }
}
```

Above is an example of how these subroutines work. 

**Key Concepts:**

1. **Fortran's `c_loc` and `c_f_pointer`:**  
   - `c_loc`: Used to obtain the C address of a Fortran array, allowing it to be passed to C/C++.
   - `c_f_pointer`: Used to convert a C pointer back into a Fortran pointer, enabling further manipulation in Fortran.

2. **C++ `extern "C"`:**
   - The `extern "C"` linkage specifies that the function uses C-style linkage, which prevents name mangling, making it callable from Fortran.

**Steps:**

1. **Fortran Side:**
   - Define a Fortran module `InteropModule` with a derived type `ArrayStruct` that represents the array's C structure.
   - Use `c_loc` to obtain a C pointer to a Fortran array and store it in the `ArrayStruct`.
   - Call the C++ subroutine to modify the array.

2. **C++ Side:**
   - Define a corresponding C++ structure `ArrayStruct` that mirrors the Fortran type.
   - Implement a function `C_subroutine` that takes a pointer to `ArrayStruct`, prints the array elements, and modifies them.
   
3. **Compilation and Linking:**
   - Compile the Fortran and C++ code, ensuring correct linkage between the two. Use `gfortran` for Fortran and `g++` for C++ with appropriate flags to ensure the Fortran program can call the C++ subroutine.

**Running the Program:**

- Upon execution, the Fortran program initializes an array, passes it to C++ for processing, and retrieves the modified array to display it. The C++ code adds 10 to each element of the array, and the changes are visible in Fortran after the subroutine call.

### Adding Shared Variables between Fortran and C++ 

#### Case 1: Single Variable

1. **C++ Side**:  
   - Add the variable declaration in `chromosphere.hpp` using the `extern` keyword.  
   - Define and use the variable in `chromosphere.cpp`. This ensures it is declared in the header for linkage but defined in the implementation file.

2. **C Side**:  
   - Add the variable to the `struct Chromo_params`.  
   - Implement inter-copy logic in both `chromo2params` and `params2chromo` to ensure the data is correctly passed between C and Fortran.

3. **Fortran Side 1 (C-Interop)**:  
   - Add the variable to `type, bind(C) :: chromo_params_C` in Fortran.  
   - The variable's type and order should exactly match the corresponding entry in `struct Chromo_params` to maintain consistency between Fortran and C.

4. **Fortran Side 2 (SWMF Format)**:  
   - Add the variable to `type :: chromo_params_F`.  
   - This follows the SWMF-specific format. The structure doesn’t need to mirror the `type, bind(C) :: chromo_params_C`, but the data must be synchronized correctly.

5. **Fortran Side 3 (Synchronization)**:  
   - Add inter-copy logic to `subroutine Chromo_SyncParamDim1Var` to ensure the new variable is synchronized between `chromo_params_F` and `chromo_params_C`.

6. **Fortran Side 4 (Initialization)**:  
   - Initialize the new variable in `subroutine Chromo_init` to set default values or prepare it for further use.

---

#### Case 2: Array Variable

Note: Only 1-D arrays (dimension(:)) are transferable between Fortran and C. Multi-dimensional arrays must be flattened to 1-D in `struct Chromo_params` and `type, bind(C) :: chromo_params_C`.

1. **C++ Side**:  
   - Declare the array as a `Vec` type in `chromosphere.hpp` using the `extern` keyword.  
   - Define and use the array in `chromosphere.cpp`.

2. **C Side**:  
   - Add a pointer to the array in `struct Chromo_params`, representing the 1-D array in C.  
   - Implement inter-copy logic in `chromo2params()` and `params2chromo()` to ensure correct data transfer.  
   - Manage memory allocation and deallocation for the array within `Chromo_allocate_params_C()` and its corresponding deallocation function.

3. **Fortran Side 1 (C-Interop)**:  
   - Add the array to `type, bind(C) :: chromo_params_C`.  
   - Use `type(c_ptr)` for array pointers in Fortran, and ensure the variable order matches `struct Chromo_params` in C.

4. **Fortran Side 2 (SWMF Format)**:  
   - Add the array to `type :: chromo_params_F` according to the SWMF format, which can differ from the C counterpart in structure but must maintain synchronized data.

5. **Fortran Side 3 (Memory Management)**:  
   - Add array allocation and deallocation in `subroutine Chromo_AllocateParams_F` and `subroutine Chromo_DeallocateParams_F`, ensuring that Fortran manages the memory efficiently.

6. **Fortran Side 4 (Data Conversion)**:  
   - Add inter-copy logic for the array in `subroutine Chromo_ConvertParam_C2F` and `subroutine Chromo_ConvertParam_F2C`.  
   - Refer to existing code for examples of how to transfer arrays between C and Fortran formats.

7. **Fortran Side 5 (Initialization)**:  
   - Initialize the array in `subroutine Chromo_init`, setting up the memory and default values as needed.
