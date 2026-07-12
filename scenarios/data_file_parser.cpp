#include "data_file_parser.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace chromosphere {

namespace {

// Trim leading/trailing ASCII whitespace in-place.
std::string trim(std::string s) {
    const auto isws = [](unsigned char c) { return c == ' ' || c == '\t' || c == '\r'; };
    while (!s.empty() && isws((unsigned char)s.back()))  s.pop_back();
    std::size_t i = 0;
    while (i < s.size() && isws((unsigned char)s[i])) ++i;
    return s.substr(i);
}

bool is_skippable(const std::string& line) {
    const auto t = trim(line);
    return t.empty() || t[0] == '#';
}

bool is_section_header(const std::string& line, std::string* name_out) {
    const auto t = trim(line);
    if (t.size() < 3 || t.front() != '[' || t.back() != ']') return false;
    if (name_out) *name_out = t.substr(1, t.size() - 2);
    return true;
}

float parse_float(const std::string& s, const std::string& ctx) {
    try {
        std::size_t pos = 0;
        const float v = std::stof(s, &pos);
        return v;
    } catch (const std::exception&) {
        throw std::runtime_error("scenario data file: bad float '" + s + "' (" + ctx + ")");
    }
}

arma::uword parse_uword(const std::string& s, const std::string& ctx) {
    try {
        std::size_t pos = 0;
        const long long v = std::stoll(s, &pos);
        if (v < 0) throw std::runtime_error("negative integer '" + s + "' (" + ctx + ")");
        return static_cast<arma::uword>(v);
    } catch (const std::exception&) {
        throw std::runtime_error("scenario data file: bad integer '" + s + "' (" + ctx + ")");
    }
}

// Read meta key=value pairs. Returns true if ns was found.
bool parse_meta_block(std::ifstream& in, ScenarioDataFile& out, std::string& next_section) {
    bool ns_seen = false;
    std::string line;
    while (std::getline(in, line)) {
        if (is_skippable(line)) continue;
        std::string sec;
        if (is_section_header(line, &sec)) { next_section = sec; return ns_seen; }
        const auto t = trim(line);
        const auto eq = t.find('=');
        if (eq == std::string::npos)
            throw std::runtime_error("scenario data file: malformed META line '" + t + "'");
        const auto k = trim(t.substr(0, eq));
        const auto v = trim(t.substr(eq + 1));
        if      (k == "ns")                 { out.ns = parse_uword(v, "META ns"); ns_seen = true; }
        else if (k == "g_si")               { out.g_si = parse_float(v, "META g_si"); }
        else if (k == "B_outer_T")          { out.B_outer_T = parse_float(v, "META B_outer_T"); }
        else if (k == "phi_g_offset_Jpkg")  { out.phi_g_offset_Jpkg = parse_float(v, "META phi_g_offset"); }
        else if (k == "topology")           { out.topology = v; }
        else if (k == "loop_half_length_m") { out.loop_half_length_m = parse_float(v, "META loop_half_length_m"); }
        // unknown keys are silently ignored — keeps the format extensible.
    }
    next_section.clear();
    return ns_seen;
}

void parse_cells_block(std::ifstream& in, ScenarioDataFile& out, std::string& next_section) {
    if (out.ns == 0)
        throw std::runtime_error("scenario data file: [CELLS] before ns= in [META]");

    out.ds_m       .set_size(out.ns);
    out.B_imh_T    .set_size(out.ns);
    out.B_iph_T    .set_size(out.ns);
    out.phi_g_imh  .set_size(out.ns);
    out.phi_g_iph  .set_size(out.ns);
    out.ne_im3     .set_size(out.ns);
    out.nn_im3     .set_size(out.ns);
    out.T_K        .set_size(out.ns);

    arma::uword rows_seen = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (is_skippable(line)) continue;
        std::string sec;
        if (is_section_header(line, &sec)) { next_section = sec; goto done; }
        std::istringstream ss(line);
        std::vector<std::string> toks;
        for (std::string t; ss >> t; ) toks.push_back(t);
        if (toks.size() != 9)
            throw std::runtime_error("scenario data file: [CELLS] expects 9 columns, got "
                                     + std::to_string(toks.size()) + " in row '" + trim(line) + "'");
        const arma::uword i = parse_uword(toks[0], "[CELLS] index");
        if (i != rows_seen)
            throw std::runtime_error("scenario data file: [CELLS] row index "
                                     + std::to_string(i) + " out of order, expected "
                                     + std::to_string(rows_seen));
        out.ds_m     (i) = parse_float(toks[1], "[CELLS] ds_m");
        out.B_imh_T  (i) = parse_float(toks[2], "[CELLS] B_imh_T");
        out.B_iph_T  (i) = parse_float(toks[3], "[CELLS] B_iph_T");
        out.phi_g_imh(i) = parse_float(toks[4], "[CELLS] phi_g_imh");
        out.phi_g_iph(i) = parse_float(toks[5], "[CELLS] phi_g_iph");
        out.ne_im3   (i) = parse_float(toks[6], "[CELLS] ne_im3");
        out.nn_im3   (i) = parse_float(toks[7], "[CELLS] nn_im3");
        out.T_K      (i) = parse_float(toks[8], "[CELLS] T_K");
        ++rows_seen;
    }
    next_section.clear();
done:
    if (rows_seen != out.ns)
        throw std::runtime_error("scenario data file: [CELLS] had "
                                 + std::to_string(rows_seen) + " rows, expected "
                                 + std::to_string(out.ns));
}

void parse_ghosts_block(std::ifstream& in, ScenarioDataFile& out, std::string& next_section) {
    out.ghost_B_T       .set_size(4);
    out.ghost_phi_g_Jpkg.set_size(4);
    out.ghost_ne_im3    .set_size(4);
    out.ghost_nn_im3    .set_size(4);
    out.ghost_T_K       .set_size(4);
    out.ghost_T_e_factor.set_size(4);

    const std::vector<std::string> tag_order = {"outer_0", "outer_1", "inner_0", "inner_1"};
    arma::uword rows_seen = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (is_skippable(line)) continue;
        std::string sec;
        if (is_section_header(line, &sec)) { next_section = sec; goto done; }
        std::istringstream ss(line);
        std::vector<std::string> toks;
        for (std::string t; ss >> t; ) toks.push_back(t);
        if (toks.size() != 7)
            throw std::runtime_error("scenario data file: [GHOSTS] expects 7 columns, got "
                                     + std::to_string(toks.size()) + " in row '" + trim(line) + "'");
        if (rows_seen >= tag_order.size())
            throw std::runtime_error("scenario data file: [GHOSTS] has more than 4 rows");
        if (toks[0] != tag_order[rows_seen])
            throw std::runtime_error("scenario data file: [GHOSTS] row "
                                     + std::to_string(rows_seen) + " has tag '" + toks[0]
                                     + "', expected '" + tag_order[rows_seen] + "'");
        out.ghost_B_T       (rows_seen) = parse_float(toks[1], "[GHOSTS] B_T");
        out.ghost_phi_g_Jpkg(rows_seen) = parse_float(toks[2], "[GHOSTS] phi_g");
        out.ghost_ne_im3    (rows_seen) = parse_float(toks[3], "[GHOSTS] ne");
        out.ghost_nn_im3    (rows_seen) = parse_float(toks[4], "[GHOSTS] nn");
        out.ghost_T_K       (rows_seen) = parse_float(toks[5], "[GHOSTS] T");
        out.ghost_T_e_factor(rows_seen) = parse_float(toks[6], "[GHOSTS] T_e_factor");
        ++rows_seen;
    }
    next_section.clear();
done:
    if (rows_seen != 4)
        throw std::runtime_error("scenario data file: [GHOSTS] had "
                                 + std::to_string(rows_seen) + " rows, expected 4");
}

} // namespace

arma::uword peek_ns_from_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("scenario data file: cannot open '" + path + "'");
    std::string line;
    bool in_meta = false;
    while (std::getline(in, line)) {
        if (is_skippable(line)) continue;
        std::string sec;
        if (is_section_header(line, &sec)) {
            in_meta = (sec == "META");
            if (!in_meta) continue;
            continue;
        }
        if (!in_meta) continue;
        const auto t = trim(line);
        const auto eq = t.find('=');
        if (eq == std::string::npos) continue;
        if (trim(t.substr(0, eq)) == "ns") return parse_uword(trim(t.substr(eq + 1)), "META ns");
    }
    throw std::runtime_error("scenario data file: ns= not found in [META] section of '" + path + "'");
}

ScenarioDataFile parse_scenario_data_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("scenario data file: cannot open '" + path + "'");

    ScenarioDataFile out;
    std::string line, current_section, next_section;

    // Scan until first section header.
    while (std::getline(in, line)) {
        if (is_skippable(line)) continue;
        if (is_section_header(line, &current_section)) break;
        throw std::runtime_error("scenario data file: data before first section header");
    }

    while (!current_section.empty()) {
        next_section.clear();
        if      (current_section == "META")   parse_meta_block  (in, out, next_section);
        else if (current_section == "CELLS")  parse_cells_block (in, out, next_section);
        else if (current_section == "GHOSTS") parse_ghosts_block(in, out, next_section);
        else throw std::runtime_error("scenario data file: unknown section [" + current_section + "]");
        current_section = next_section;
    }

    if (out.ns == 0)                throw std::runtime_error("scenario data file: missing [META] ns=");
    if (out.ds_m.n_elem == 0)       throw std::runtime_error("scenario data file: missing [CELLS]");
    if (out.ghost_B_T.n_elem != 4)  throw std::runtime_error("scenario data file: missing [GHOSTS]");
    return out;
}

} // namespace chromosphere
