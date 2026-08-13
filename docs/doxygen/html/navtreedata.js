/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "Chromosphere2026", "index.html", [
    [ "Two solvers, one mesh", "index.html#autotoc_md35", null ],
    [ "Where to start", "index.html#autotoc_md36", null ],
    [ "Module index", "index.html#autotoc_md37", null ],
    [ "The release configuration in one paragraph", "index.html#autotoc_md38", null ],
    [ "Primary references", "index.html#autotoc_md39", null ],
    [ "Source tree and architecture", "architecture.html", [
      [ "The dependency rule", "architecture.html#autotoc_md0", null ],
      [ "Annotated file tree", "architecture.html#autotoc_md1", [
        [ "Shared headers (repository root)", "architecture.html#autotoc_md2", null ],
        [ "Shared implementation, <span class=\"tt\">src/*.cpp</span>", "architecture.html#autotoc_md3", null ],
        [ "Release solver, <span class=\"tt\">src/single_fluid/</span>", "architecture.html#autotoc_md4", null ],
        [ "Historical two-fluid solver, <span class=\"tt\">src/two_fluid/</span>", "architecture.html#autotoc_md5", null ],
        [ "Scenarios, <span class=\"tt\">scenarios/</span>", "architecture.html#autotoc_md6", null ],
        [ "Everything else", "architecture.html#autotoc_md7", null ]
      ] ],
      [ "State ownership", "architecture.html#autotoc_md8", null ],
      [ "Solver selection", "architecture.html#autotoc_md9", null ],
      [ "Lifecycle of one run", "architecture.html#autotoc_md10", null ],
      [ "Build layout", "architecture.html#autotoc_md11", [
        [ "Registered CTest cases", "architecture.html#autotoc_md12", null ],
        [ "OpenMP", "architecture.html#autotoc_md13", null ]
      ] ]
    ] ],
    [ "Configuration reference", "configuration.html", [
      [ "Command line", "configuration.html#autotoc_md14", null ],
      [ "Environment variables", "configuration.html#autotoc_md15", [
        [ "Run control and output cadence", "configuration.html#autotoc_md16", null ],
        [ "Diagnostics and profiling", "configuration.html#autotoc_md17", null ],
        [ "Solver selection", "configuration.html#autotoc_md18", null ],
        [ "<span class=\"tt\">model_column</span> / <span class=\"tt\">model_gentle</span> model knobs (the <span class=\"tt\">ISO_*</span> family)", "configuration.html#autotoc_md19", null ],
        [ "Mesh refinement", "configuration.html#autotoc_md20", null ],
        [ "Flare, PFSS and C7 knobs", "configuration.html#autotoc_md21", null ],
        [ "OpenMP runtime", "configuration.html#autotoc_md22", null ]
      ] ],
      [ "Reference-only overrides — never use in production", "configuration.html#autotoc_md23", null ],
      [ "Rejected in release mode", "configuration.html#autotoc_md24", null ],
      [ "Validated production defaults", "configuration.html#autotoc_md25", null ]
    ] ],
    [ "Output file formats", "io_formats.html", [
      [ "Main snapshot file", "io_formats.html#autotoc_md26", [
        [ "<span class=\"tt\">state_rows == 3</span> — release snapshot", "io_formats.html#autotoc_md27", null ],
        [ "<span class=\"tt\">state_rows == 7</span> — two-fluid snapshot", "io_formats.html#autotoc_md28", null ]
      ] ],
      [ "<span class=\"tt\">&lt;out&gt;.gamma_diag</span>", "io_formats.html#autotoc_md29", null ],
      [ "<span class=\"tt\">&lt;out&gt;.faceflux</span>", "io_formats.html#autotoc_md30", null ],
      [ "<span class=\"tt\">&lt;out&gt;.outercond</span>", "io_formats.html#autotoc_md31", [
        [ "Sign and geometry conventions", "io_formats.html#autotoc_md32", null ]
      ] ],
      [ "<span class=\"tt\">outputs/output.log</span>", "io_formats.html#autotoc_md33", null ],
      [ "Reading these in Python", "io_formats.html#autotoc_md34", null ]
    ] ],
    [ "Naming conventions", "naming.html", [
      [ "Grid-location suffix", "naming.html#autotoc_md40", null ],
      [ "Variable bodies", "naming.html#autotoc_md41", null ],
      [ "Index constants", "naming.html#autotoc_md42", null ],
      [ "Function names", "naming.html#autotoc_md43", null ],
      [ "Units", "naming.html#autotoc_md44", null ]
    ] ],
    [ "Numerical method", "numerics.html", [
      [ "The release scheme", "numerics.html#autotoc_md45", null ],
      [ "Finite-volume discretization", "numerics.html#autotoc_md46", null ],
      [ "Reconstruction", "numerics.html#autotoc_md47", null ],
      [ "Predictor and numerical flux", "numerics.html#autotoc_md48", null ],
      [ "Reference-free hydrostatic balance", "numerics.html#autotoc_md49", null ],
      [ "Implicit conduction", "numerics.html#autotoc_md50", null ],
      [ "EOS inversion performance", "numerics.html#autotoc_md51", null ],
      [ "Timestep control", "numerics.html#autotoc_md52", null ],
      [ "The two-fluid scheme", "numerics.html#autotoc_md53", null ],
      [ "Rejected, and not to be reintroduced without new evidence", "numerics.html#autotoc_md54", null ],
      [ "See also", "numerics.html#autotoc_md55", null ]
    ] ],
    [ "Physical model", "physics_model.html", [
      [ "Field-aligned reduction", "physics_model.html#autotoc_md56", null ],
      [ "The release model", "physics_model.html#release_physics", [
        [ "Governing equations", "physics_model.html#autotoc_md58", null ],
        [ "Closure", "physics_model.html#autotoc_md59", null ],
        [ "Conduction", "physics_model.html#autotoc_md60", null ],
        [ "Boundaries", "physics_model.html#autotoc_md61", null ],
        [ "What the release deliberately does not contain", "physics_model.html#autotoc_md62", null ]
      ] ],
      [ "The historical two-fluid model", "physics_model.html#two_fluid_physics", [
        [ "Conserved state", "physics_model.html#autotoc_md64", null ],
        [ "Source physics", "physics_model.html#autotoc_md65", [
          [ "Ionization network", "physics_model.html#autotoc_md66", null ],
          [ "Radiative losses", "physics_model.html#autotoc_md67", null ],
          [ "Evaporation drivers", "physics_model.html#autotoc_md68", null ]
        ] ],
        [ "Known issue", "physics_model.html#autotoc_md69", null ],
        [ "See also", "physics_model.html#autotoc_md70", null ]
      ] ]
    ] ],
    [ "Scenario reference", "scenario_reference.html", [
      [ "Comparison", "scenario_reference.html#autotoc_md71", null ],
      [ "<span class=\"tt\">model_column</span> — the release scenario", "scenario_reference.html#autotoc_md72", null ],
      [ "<span class=\"tt\">model_gentle</span> — the historical two-fluid column", "scenario_reference.html#autotoc_md73", null ],
      [ "<span class=\"tt\">model_c7</span> — quiet-Sun baseline and shared C7 library", "scenario_reference.html#autotoc_md74", null ],
      [ "<span class=\"tt\">model_flare</span> — explosive evaporation", "scenario_reference.html#autotoc_md75", null ],
      [ "<span class=\"tt\">analytic_canopy</span> — exponential magnetic canopy", "scenario_reference.html#autotoc_md76", null ],
      [ "<span class=\"tt\">pfss_field_line</span> — traced potential-field line", "scenario_reference.html#autotoc_md77", null ]
    ] ],
    [ "Validation and known limitations", "validation.html", [
      [ "Running the checks", "validation.html#autotoc_md78", null ],
      [ "What the C++ suite covers", "validation.html#autotoc_md79", null ],
      [ "Registered CTest cases", "validation.html#autotoc_md80", null ],
      [ "Known limitations", "validation.html#autotoc_md81", null ],
      [ "Rejected approaches", "validation.html#autotoc_md82", null ],
      [ "Evidence documents", "validation.html#autotoc_md83", null ]
    ] ],
    [ "Topics", "topics.html", "topics" ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Functions", "namespacemembers_func.html", null ],
        [ "Variables", "namespacemembers_vars.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", null ],
        [ "Variables", "functions_vars.html", "functions_vars" ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"analytic__canopy_8cpp.html",
"group__runtime.html#gga158464c6b9ead9c4a6c5b5e5dd586ba3af2c5bf353075e0b4d3b6bbe5c0d061be",
"structchromosphere_1_1_grid.html#a173b0a898d8bbe099b59ace0072ef850",
"structchromosphere_1_1_scenario_data_file.html"
];

const SYNCONMSG = 'click to disable panel synchronization';
const SYNCOFFMSG = 'click to enable panel synchronization';
const LISTOFALLMEMBERS = 'List of all members';