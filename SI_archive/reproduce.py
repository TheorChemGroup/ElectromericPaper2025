import re
import csv
import subprocess
import sys
import os
import math
from itertools import zip_longest

class GoodVibesManager:

    '''Verifies the installation of GoodVibes2 from https://github.com/TheorChemGroup/GoodVibes2.git'''

    def __init__(self):
        pass

    def compare_version(self, line = "GoodVibes-2 v1.0") -> bool:
        '''Compares version of goodvibes called by 'python -m goodvibes' with GoodVibes2 version'''
        try:
            result = subprocess.run(
                ["python", "-m", "goodvibes", "-h"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            
            output = result.stdout
            match = re.search(r'--version\s+([^\s]+)', output)
            
            if match is None:
                print("Could not find version info in goodvibes output.")

            elif match.group(1) == line:
                print("✓ Versions match.")
                return True
            
            else:
                print("✗ Versions do not match. Make sure to install Goodvibes2 from https://github.com/TheorChemGroup/GoodVibes2.git\n")
                return False

        except FileNotFoundError:
            print("Error: 'goodvibes' command not found. Is GoodVibes2 installed and in your PATH?")
            return False
        except subprocess.CalledProcessError as e:
            print(f"Error running 'goodvibes --version': {e}")
            print(f"stderr: {e.stderr}")
            return False


class GoodVibesRunner:

    '''Runs GoodVibes2 to compute thermochemical data from electronic structure calculations in ORCA.'''

    def __init__(self):
        pass

    def run(self, filepath, temperature=298.15): 
        '''Runs GoodVibes with the following settings: 
        Truhlar approximation for the quasy-harmonic entropy correction,
        frequency entropic cut-off of 175 cm^-1, temperature of 298.15 K, 
        concentration of 1 mol/l, scaling factor for vibrational frequencies of 1, 
        conversion of imaginary frequencies into real ones'''
        cmd = [
            "python",
            "-m",
            "goodvibes",
            filepath,
            "--qs", "truhlar",
            "--fs", "175",
            "-t", str(temperature),
            "--no_plot",
            "-c", "1",
            "-v", "1.0",
            "--invertifreq", "auto",
        ]
        try:
            env = os.environ.copy()
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                env=env,
            )
            if result.returncode != 0:
                print(f"GoodVibes error for {filepath}:\n{result.stderr}")
                sys.exit()
                return None, None
            

            return self.extract_energies_from_goodvibes_output(result.stdout)

        except Exception as e:
            print(f"Error running GoodVibes on {filepath}: {e}")
            sys.exit()
            return None, None
        
    def extract_energies_from_goodvibes_output(self, output):
        '''Extracts Electronic (E) and Gibbs energies with quasi-harmonic corrections (qhG) 
        from GoodVibes2 output lines'''
        lines = output.splitlines()
        header_line = None
        data_line = None

        for line in lines:
            if "E        ZPE             H        T.S     T.qh-S          G(T)          qh-G(T)" in line:
                header_line = line
            elif line.strip().startswith('o'):
                data_line = line
        try:
            headers = header_line.split()
            values = data_line.split()
            E_index = headers.index('E') + 1
            qhg_index = headers.index('qh-G(T)') + 1
            E_value = float(values[E_index])
            qhg_value = float(values[qhg_index])

            return E_value, qhg_value
        
        except (ValueError, IndexError) as e:
            print(f"Warning: Error parsing: {e}")
            return None, None


class LogFileParser:

    '''Obtains total charge, multiplicity and number of imaginary frequencies from ORCA .log files'''

    TOTAL_CHARGE_LINE = "Total Charge"
    MULTIPLICITY_LINE = "Multiplicity"

    def __init__(self, filepath):
        self.filepath = filepath
        self.charge = None
        self.multiplicity = None
        self.num_imfreq = 0
    
    def parse(self):
        '''Returns charge, multiplicity and number of imaginary frequencies obtained from ORCA .log file'''
        try:
            with open(self.filepath, 'r', encoding='utf-8') as file:
                lines = file.readlines()
                self.parse_charge_and_multiplicity(lines)
                self.parse_number_of_imaginary_modes(lines)
        except Exception as e:
            print(f"error reading {self.filepath}: {e}")
        return self.charge, self.multiplicity, self.num_imfreq
    
    def parse_charge_and_multiplicity(self, lines):
        '''Returns charge, multiplicity obtained from ORCA .log file'''
        for line in lines:
            if self.TOTAL_CHARGE_LINE in line:
                parts = line.strip().split()
                for part in reversed(parts):
                    if part.lstrip('-').isdigit():
                        self.charge = int(part)
                        break
            if self.MULTIPLICITY_LINE in line:
                parts = line.strip().split()
                for part in reversed(parts):
                    if part.isdigit():
                        self.multiplicity = int(part)
                        break
            if self.charge is not None and self.multiplicity is not None:
                break

    def parse_number_of_imaginary_modes(self, lines):
        '''Returns charge number of imaginary frequencies obtained from ORCA .log file'''
        vib_sections = []
        current_section = None
        for i, line in enumerate(lines):
            if "VIBRATIONAL FREQUENCIES" in line:
                current_section = {"start": i, "end": None, "count": 0}
            elif current_section and "NORMAL MODES" in line:
                current_section["end"] = i
                vib_sections.append(current_section)
                current_section = None
            elif current_section and "***imaginary mode***" in line:
                current_section["count"] += 1
        if vib_sections:
            self.num_imfreq = vib_sections[-1]["count"]


class Utilities:

    @staticmethod
    def get_dir(dir_name, root_name):
        '''Returns absolute path of directory inside root directory'''
        root = Utilities.find_root_dir(root_name)
        if isinstance(dir_name, str):
            path = os.path.join(root, dir_name)
        else:
            path = os.path.join(root, *dir_name)
        if not os.path.isdir(path):
            raise FileNotFoundError(f"Directory '{path}' does not exist. Make sure to call the script from SI_archive folder")
        return path
    
    @staticmethod
    def find_root_dir(root_name):
        '''Returns absolute path of root directory'''
        current_dir = os.path.abspath(os.getcwd())
        while True:
            if os.path.basename(current_dir)==root_name:
                return current_dir
            parent_dir = os.path.dirname(current_dir)
            if parent_dir==current_dir:
                raise FileNotFoundError(f"Root directory '{root_name}' not found. Make sure to call the script from SI_archive folder")
            current_dir = parent_dir

    @staticmethod
    def extract_substituent_from_filename(filename, subdir):
        '''Extracts substituent name from filename'''
        base = os.path.basename(filename)
        if base.endswith('.log'):
            base = base[:-4]
        
        if base.endswith("PhenO") or base.endswith("PhenOH"):
            return base.split("_")[0]
        
        if subdir and base.startswith(subdir + "_"):
            return base[len(subdir) + 1:]
        
        return base

    @staticmethod
    def write_summary_csv(data, filename):
        '''Basic function for writing data into .csv file'''
        if not data:
            return
        fieldnames = list(data[0].keys())
        with open(filename, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
    

class EnergyProfileAnalyzer:

    '''Calculates energies of states with respect to 3'.
    Saves table with all general structure information to .csv'''

    ROOT_DIR_NAME = "SI_archive"
    ENERGY_PROFILES_DIR = ["Energy_profiles"]
    OUTPUT_CSV = "energy_profiles.csv"

    def __init__(self):
        self.root_dir = Utilities.find_root_dir(self.ROOT_DIR_NAME)
        self.energy_profiles_dir = Utilities.get_dir(self.ENERGY_PROFILES_DIR, self.ROOT_DIR_NAME)
        self.goodvibes_manager = GoodVibesManager()
        self.goodvibes_runner = GoodVibesRunner()

    def extract_all_data(self, dir_name, root_name, goodvibes_runner):
        '''Extracts information about file, charge, multiplicity, number of imaginary frequencies,
        Electronic energy and Gibbs energy with quasi-harmonic corrections for all files in Energy_profiles directory.'''
        base_path = Utilities.get_dir(dir_name, root_name)
        all_files_data = []

        for subdir_name in sorted(os.listdir(base_path)):
            subdir_path = os.path.join(base_path, subdir_name)
            if not os.path.isdir(subdir_path):
                continue
                
            for filename in sorted(os.listdir(subdir_path)):
                if not filename.endswith(".log"):
                    continue
                filepath = os.path.join(subdir_path, filename)
                parser = LogFileParser(filepath)
                charge, mult, num_imfreq = parser.parse()
                substituent = Utilities.extract_substituent_from_filename(filename, subdir_name)
                print(f"Processing file {filename}")  
                
                E, qhG = goodvibes_runner.run(filepath)

                entry = {
                    "Name": filename,
                    "Subdir": subdir_name,
                    "Substituent": substituent,
                    "Charge": charge,
                    "Multiplicity": mult,
                    "NumImFreq": num_imfreq,
                    "Energy": E,
                    "qhG": qhG,
                    "FilePath": filepath
                }
                all_files_data.append(entry)

        return all_files_data
    
    def calculate_relative_energies(self, all_data):
        '''Calculates Electronic and qhG energies relative to 3'.'''
        reference_energy = {}
        reference_qhG = {}
        pyridine_energy = None
        pyridine_qhG = None

        for entry in all_data:
            if entry.get("Subdir") == "4'" and entry.get("Substituent") == "pyridine":
                pyridine_energy = entry.get("Energy")
                pyridine_qhG = entry.get("qhG")
                break

        for entry in all_data:
            if entry.get("Subdir") == "3'":
                sub = entry.get("Substituent")
                if sub:
                    if entry.get("Energy") is not None:
                        reference_energy[sub] = entry["Energy"]
                    if entry.get("qhG") is not None:
                        reference_qhG[sub] = entry["qhG"]

        results = []
        for entry in all_data:
            sub = entry.get("Substituent")
            if not sub or sub == "pyridine":
                continue
            current_energy = entry.get("Energy")
            current_qhG = entry.get("qhG")
            subdir = entry.get("Subdir")
            if subdir == "4'":
                if pyridine_energy is None or pyridine_qhG is None:
                    continue
                energy = current_energy + pyridine_energy if current_energy is not None else None
                qhG = current_qhG + pyridine_qhG if current_qhG is not None else None
            else:
                energy = current_energy
                qhG = current_qhG
            e_rel = (energy - reference_energy[sub]) if (energy is not None and sub in reference_energy) else None
            qhG_rel = (qhG - reference_qhG[sub]) if (qhG is not None and sub in reference_qhG) else None
            e_rel_kcal = e_rel * 627.5095 if e_rel is not None else None
            qhG_rel_kcal = qhG_rel * 627.5095 if qhG_rel is not None else None

            results.append({
                "Name": entry["Name"],
                "Subdir": entry["Subdir"],
                "Substituent": sub,
                "Charge": entry.get("Charge", ""),
                "Multiplicity": entry.get("Multiplicity", ""),
                "NumImFreq": entry.get("NumImFreq", 0),
                "Energy (Eh)": f"{energy:.5f}" if energy is not None else "",
                "qhG (Eh)": f"{qhG:.5f}" if qhG is not None else "",
                "E rel (kcal/mol)": f"{e_rel_kcal:.5f}" if e_rel_kcal is not None else "",
                "qhG rel (kcal/mol)": f"{qhG_rel_kcal:.5f}" if qhG_rel_kcal is not None else "" 
            })
        return results

    def get_energy_profiles(self, data):
        '''Creates Fig_S20.csv and Fig_2E.csv based on calculated relative qhG'''
        substituents = set()
        group_A = ["Me", "4-OMePh", "4-NO2Ph", "NO2"]
        group_B = ["H", "CHCH2", "CF3", "CN"]
        group_C = ["F", "N3", "CCH"]
        group_D = ["NMe2", "CHCHNMe2", "CO2Me", "CHCHCO2Me"]
        for entry in data:
            if entry.get("Substituent"):
                substituents.add(entry["Substituent"])
        if not substituents:
            print("Warning: No substituents found in data")
            return
        plot_data = {}
        for entry in data:
            sub = entry.get("Substituent")
            subdir = entry.get("Subdir")
            qhG_rel_kcal = entry.get("qhG rel (kcal/mol)")
            
            if sub and subdir and qhG_rel_kcal and qhG_rel_kcal.strip():
                try:
                    energy = float(qhG_rel_kcal)
                    if sub not in plot_data:
                        plot_data[sub] = {}
                    plot_data[sub][subdir] = energy
                except ValueError:
                    continue
        
        PROFILE_ORDER = ["3'", "TS1", "IM", "TS2", "4'"]
        
        substituents = sorted(plot_data.keys())

        with open("Fig_S20A.csv", 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Substituent"] + PROFILE_ORDER)
            for sub in substituents:
                if sub in group_A:
                    row = [sub]
                    for stage in PROFILE_ORDER:
                        energy = plot_data[sub].get(stage)
                        row.append(f"{energy:.2f}" if energy is not None else '')
                    writer.writerow(row)

        with open("Fig_S20B.csv", 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Substituent"] + PROFILE_ORDER)
            for sub in substituents:
                if sub in group_B:
                    row = [sub]
                    for stage in PROFILE_ORDER:
                        energy = plot_data[sub].get(stage)
                        row.append(f"{energy:.2f}" if energy is not None else '')
                    writer.writerow(row)

        with open("Fig_S20C.csv", 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Substituent"] + PROFILE_ORDER)
            for sub in substituents:
                if sub in group_C:
                    row = [sub]
                    for stage in PROFILE_ORDER:
                        energy = plot_data[sub].get(stage)
                        row.append(f"{energy:.2f}" if energy is not None else '')
                    writer.writerow(row)

        with open("Fig_S20D.csv", 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Substituent"] + PROFILE_ORDER)
            for sub in substituents:
                if sub in group_D:
                    row = [sub]
                    for stage in PROFILE_ORDER:
                        energy = plot_data[sub].get(stage)
                        row.append(f"{energy:.2f}" if energy is not None else '')
                    writer.writerow(row)

        fig_2E = "Fig_2E.csv"
        with open(fig_2E, 'w', newline='') as fig_2C_file:
            writer = csv.writer(fig_2C_file)
            writer.writerow(["Substituent"] + PROFILE_ORDER)
            for sub in substituents:
                if sub in ["NO2", "4-NO2Ph"]:
                    row = [sub]
                    for stage in PROFILE_ORDER:
                        energy = plot_data[sub].get(stage)
                        row.append(f"{energy:.2f}" if energy is not None else '')
                    writer.writerow(row)            
        print(f"\nCreated {fig_2E}, Fig_S20A.csv, Fig_S20B.csv, Fig_S20C.csv, Fig_S20D.csv\n")

    def run_analysis(self):
        '''Runs analysis to create .csv files from Energy profiles data'''
        self.goodvibes_manager.compare_version()
        all_data = self.extract_all_data(self.ENERGY_PROFILES_DIR, EnergyProfileAnalyzer.ROOT_DIR_NAME, self.goodvibes_runner)
        energy_profile_data = self.calculate_relative_energies(all_data)
        self.get_energy_profiles(energy_profile_data)
        Utilities.write_summary_csv(energy_profile_data, "Structural_summary_data.csv")
        return energy_profile_data


class OrbitalEnergyExtractor:

    '''Extracts orbital energies used for E_HU calculations,
    Parses basic info of simplified_IM structures to orbital_data.csv,
    Saves tables for correlations Fig.3E and Fig3G to .csv .'''

    ROOT_DIR_NAME = "SI_archive"
    E_HU_DIR = ["E_HU"]

    def __init__(self):
        self.root_dir = Utilities.find_root_dir(self.ROOT_DIR_NAME)
        self.e_hu_dir = os.path.join(self.root_dir, *self.E_HU_DIR)

    def extract_orbital_data(self):
        '''Extracts file information, charge, multiplicity, HOMO, LUMO and LUMO+1 orbital energies
        and HOMO-LUMO (eV) Gap from files in E_HU directory. Writes data to Fog_S20.csv.'''
        data = []
        for filename in sorted(os.listdir(self.e_hu_dir)):
            if not (filename.startswith("simplified_IM_") and filename.endswith(".log")):
                continue
            filepath = os.path.join(self.e_hu_dir, filename)
            substituent = filename.replace("simplified_IM_", "").replace(".log", "")
            print(f"Processing {filename}", end=" \n")
            try:
                logfile_parser = LogFileParser(filepath)
                charge, mult, num_imfreq = logfile_parser.parse()
                homo, lumo, lumo_plus_1, gap = self.parse_orbital_energies(filepath)
                if None in [homo, lumo]:
                    print("Failed - missing orbital data")
                    continue
                data.append({
                    "Filename": filename,
                    "Substituent": substituent,
                    "Charge": charge,
                    "Multiplicity": mult,
                    "NumImFreq": num_imfreq,
                    "HOMO (eV)": homo,
                    "LUMO (eV)": lumo,
                    "LUMO+1 (eV)": lumo_plus_1,
                    "HOMO-LUMO Gap (eV)": gap
                })
            except Exception as e:
                print(f"Error: {str(e)}")
                continue
        if data:
            orbital_data = "Fig_S21.csv"
            csv_file_path = os.path.join("", orbital_data)
            with open(csv_file_path, "w", newline="", encoding="utf-8") as csvfile:
                fieldnames = ["Filename", "Substituent", "Charge", "Multiplicity", "NumImFreq",
                            "HOMO (eV)", "LUMO (eV)", "LUMO+1 (eV)", "HOMO-LUMO Gap (eV)"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for row in data:
                    writer.writerow(row)
            print(f"\nCreated {orbital_data}\n")
        else:
            print("No orbital data extracted, CSV file not created.")
        return data

    def parse_orbital_energies(self, filepath):
        '''Returns HOMO, LUMO and LUMO+1 orbital energies
        and calculates HOMO-LUMO (eV) Gap from a .log file in E_HU directory.'''
        with open(filepath, 'r', encoding='utf-8') as file:
            lines = file.readlines()[-200000:] 
        orbital_section_start = None
        for i in range(len(lines) - 1, -1, -1):
            if "ORBITAL ENERGIES" in lines[i]:
                orbital_section_start = i
                break
        if orbital_section_start is None:
            return None, None, None, None 

        homo = None
        lumo = None
        lumo_plus_1 = None
        found_lumo = False
        unoccupied_count = 0

        for line in lines[orbital_section_start:]:
            if "Alpha spin orbitals" in line:
                break
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                occupation = float(parts[1])
                energy = float(parts[3])
                if occupation > 0.5:
                    homo = energy
                elif not found_lumo:
                    lumo = energy
                    found_lumo = True
                    unoccupied_count += 1
                else:
                    unoccupied_count += 1
                    if unoccupied_count == 2:
                        lumo_plus_1 = energy
                        break
            except (ValueError, IndexError):
                continue

        gap = lumo - homo if (homo is not None and lumo is not None) else None
        return homo, lumo, lumo_plus_1, gap


    def save_to_csv(self, orbital_data, energy_profile_data):
        """Saves correlations of E_HU with G_IM and G_TS2 to .csv files.
        E_HU is equal to HOMO-LUMO Gap in most cases (see Fig_S21 in Supplementary Information). 
        In two cases (R=N3 and R=4-NO2Ph), E_HU is equal to HOMO-LUMO+1 Gap.
        Groups substituents based on their Donating/Accepting ability and High or Low polarizability and
        creates Fig_3A.csv, Fig_3B.csv, Fig_3C.csv, Fig_3E.csv and Fig_3G.csv"""
        if not orbital_data:
            print("No data to save")
            return

        fig3E_data = []
        fig3G_data = []
        fig3A_data = []

        PI_CONJUGATED_SUBS = {
            "CHCH2", "CHCHCO2Me", "4-NO2Ph", "CHCHNMe2",
            "CO2Me", "CCH", "4-OMePh", "NO2", "N3", "CN"
        }
        DONORS = {"CHCHNMe2", "CHCH2", "4-OMePh", "N3",
                "Me", "NMe2", "H", "F"
        }

        ts2_groups = {
            "donors": [],
            "acceptors": [],
            "high (pi-conj)": [],
            "low (non-pi-conj)": []
        }
        im_groups = {
            "donors": [],
            "acceptors": [],
            "high (pi-conj)": [],
            "low (non-pi-conj)": []
        }

        for entry in orbital_data:
            sub = entry['Substituent'].strip()  
            
            ts2_entry = next((e for e in energy_profile_data 
                            if e['Substituent'].strip() == sub and e['Subdir'] == "TS2"), None)
            im_entry = next((e for e in energy_profile_data 
                            if e['Substituent'].strip() == sub and e['Subdir'] == "IM"), None)
            
            try:
                e_im = float(im_entry['qhG rel (kcal/mol)']) if im_entry else None
            except (TypeError, ValueError, KeyError):
                e_im = None

            try:
                ts2_energy = float(ts2_entry['qhG rel (kcal/mol)']) if ts2_entry else None
            except (TypeError, ValueError):
                ts2_energy = None

            if sub in ["N3", "4-NO2Ph"]:
                lumo_plus_1 = entry.get("LUMO+1 (eV)", 0)
                homo = entry.get("HOMO (eV)", 0)
                if isinstance(lumo_plus_1, (int, float)) and isinstance(homo, (int, float)):
                    orbital_energy = lumo_plus_1 - homo
                else:
                    orbital_energy = None
            else:
                orbital_energy = entry.get("HOMO-LUMO Gap (eV)")

            orbital_energy_fmt = f"{orbital_energy:.2f}" if isinstance(orbital_energy, (int, float)) else ""

            e_im_fmt = f"{e_im:.2f}" if isinstance(e_im, (int, float)) else ""
            ts2_energy_fmt = f"{ts2_energy:.2f}" if isinstance(ts2_energy, (int, float)) else ""

            fig3E_data.append({
                "Substituent": sub,
                "E_HU (eV)": orbital_energy_fmt,
                "G_IM (kcal/mol)": e_im_fmt
            })

            fig3G_data.append({
                "Substituent": sub,
                "E_HU (eV)": orbital_energy_fmt,
                "G_TS2 (kcal/mol)": ts2_energy_fmt
            })

            group_high_low = "high (pi-conj)" if sub in PI_CONJUGATED_SUBS else "low (non-pi-conj)"
            group_d_a = "donors" if sub in DONORS else "acceptors"

            fig3A_data.append({
                "Substituent": sub,
                "Group": group_high_low,
                "G_IM (kcal/mol)": e_im_fmt,
                "G_TS2 (kcal/mol)": ts2_energy_fmt
            })

            if ts2_energy is not None:
                ts2_groups[group_d_a].append(f"{ts2_energy:.2f}")
                ts2_groups[group_high_low].append(f"{ts2_energy:.2f}")
            if e_im is not None:
                im_groups[group_d_a].append(f"{e_im:.2f}")
                im_groups[group_high_low].append(f"{e_im:.2f}")

        with open("Fig_3C.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            headers = ["donors", "acceptors", "high (pi-conj)", "low (non-pi-conj)"]
            writer.writerow(headers)
            for row in zip_longest(ts2_groups["donors"], ts2_groups["acceptors"],
                                ts2_groups["high (pi-conj)"], ts2_groups["low (non-pi-conj)"],
                                fillvalue=""):
                writer.writerow(row)
        print("\nCreated Fig_3C.csv\n")

        with open("Fig_3B.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            for row in zip_longest(im_groups["donors"], im_groups["acceptors"],
                                im_groups["high (pi-conj)"], im_groups["low (non-pi-conj)"],
                                fillvalue=""):
                writer.writerow(row)
        print("\nCreated Fig_3B.csv\n")

        with open("Fig_3E.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Substituent", "E_HU (eV)", "G_IM (kcal/mol)"])
            for row in fig3E_data:
                writer.writerow([row["Substituent"], row["E_HU (eV)"], row["G_IM (kcal/mol)"]])
        print("\nCreated Fig_3E.csv\n")

        with open("Fig_3G.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Substituent", "E_HU (eV)", "G_TS2 (kcal/mol)"])
            for row in fig3G_data:
                writer.writerow([row["Substituent"], row["E_HU (eV)"], row["G_TS2 (kcal/mol)"]])
        print("\nCreated Fig_3G.csv\n")

        with open("Fig_3A.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Substituent", "Group", "G_IM (kcal/mol)", "G_TS2 (kcal/mol)"])
            for row in fig3A_data:
                writer.writerow([
                    row["Substituent"],
                    row["Group"],
                    row["G_IM (kcal/mol)"],
                    row["G_TS2 (kcal/mol)"]
                ])
        print("\nCreated Fig_3A.csv\n")


class HammettConstants:
    '''Obtains values of sigma_p_minus for several substituents based on sigma_p_minus formula (see Computational details in SI)'''
    ROOT_DIR_NAME = "SI_archive"
    HAMMETT_DIR = "Hammett_constants"
    R = 1.9872036e-3 
    T = 298.15       
    LN10 = math.log(10)

    def __init__(self, energy_profiles_data):
        self.energy_profiles_data = energy_profiles_data
        self.root_dir = Utilities.find_root_dir(self.ROOT_DIR_NAME)
        self.hammett_dir = os.path.join(self.root_dir, self.HAMMETT_DIR)
        self.goodvibes_manager = GoodVibesManager()
        self.goodvibes_runner = GoodVibesRunner()

    def extract_sub_and_suffix(self, filename):
        """
        Extracts state of the molecule (PhenO- or PhenOH) and substituent in its para-position.
        """
        base = os.path.basename(filename)
        if base.endswith('.log'):
            base = base[:-4]
        m = re.match(r"(?:(.+)_)?(PhenO|PhenOH)$", base)
        if m:
            sub = m.group(1) if m.group(1) else ""  
            suffix = m.group(2)
            return sub, suffix
        return None, None

    def parse_hammett_data(self):
        """
        Parses all log files in Hammett_constants directory.
        Returns a list of entries with general data (Subdir is always "").
        """
        entries = []
        for filename in sorted(os.listdir(self.hammett_dir)):
            if not filename.endswith(".log"):
                continue
                
            filepath = os.path.join(self.hammett_dir, filename)
            
            substituent, sub = self.extract_sub_and_suffix(filename)
            
            parser = LogFileParser(filepath)
            charge, mult, num_imfreq = parser.parse()
            
            print(f"Processing file {filename}")  
            E, qhG = self.goodvibes_runner.run(filepath)
            
            if E is None or qhG is None:
                print(f"Warning: GoodVibes failed for {filename}")
                continue

            entry = {
                "Name": filename,
                "Substituent": substituent,
                "Charge": charge,
                "Multiplicity": mult,
                "NumImFreq": num_imfreq,
                "Energy": E,
                "qhG": qhG,
                "FilePath": filepath
            }
            entries.append(entry)
        
        return entries

    def save_hammett_data_to_csv(self, hammett_entries):
        """Saves the parsed Hammett constants data to a CSV file"""
        if not hammett_entries:
            print("No Hammett data to save")
            return
    
        csv_data = []
        for entry in hammett_entries:
            csv_data.append({
                "Filename": entry["Name"],
                "Substituent": entry["Substituent"],
                "Charge": entry["Charge"],
                "Multiplicity": entry["Multiplicity"],
                "NumImFreq": entry["NumImFreq"],
                "Energy (Eh)": f"{entry['Energy']:.8f}" if entry['Energy'] is not None else "",
                "qhG (Eh)": f"{entry['qhG']:.8f}" if entry['qhG'] is not None else ""
            })
        
        Utilities.write_summary_csv(csv_data, "Hammett_data_summary.csv")
        print("\nCreated Hammett_data_summary.csv\n")

    def calculate_sigma_p(self, hammett_entries):
        """
        Calculates sigma_p_minus for each substituent X using:
        σ_X = (qhG_(XPhenOH) - qhG_(XPhenO) - (qhG_PhenOH - qhG_PhenO)) / (R * T * ln(10))
        """
        PhenO = "PhenO"
        PhenOH = "PhenOH"
        qhG_lookup = {}
        for entry in hammett_entries:
            sub, suffix = self.extract_sub_and_suffix(entry["Name"])
            if sub is not None and suffix is not None:
                qhG_lookup[(sub, suffix)] = entry["qhG"]*627.5095

        qhG_PhOH = qhG_lookup.get(("", PhenOH))
        qhG_PhO = qhG_lookup.get(("", PhenO))
        if qhG_PhOH is None or qhG_PhO is None:
            raise ValueError("Reference PhenOH or PhenO data missing in Hammett constants folder.")

        sigma_data = []
        all_subs = set(sub for (sub, suffix) in qhG_lookup.keys() if sub and suffix in (PhenOH, PhenO))
        for X in sorted(all_subs):
            qhG_XPhOH = qhG_lookup.get((X, PhenOH))
            qhG_XPhO = qhG_lookup.get((X, PhenO))
            if qhG_XPhOH is None or qhG_XPhO is None:
                continue  

            numerator = (qhG_XPhOH - qhG_XPhO) - (qhG_PhOH - qhG_PhO)
            denominator = self.R * self.T * self.LN10
            sigma_p_minus = numerator / denominator
            sigma_data.append({
                "Substituent": X,
                "sigma_p_minus": sigma_p_minus
            })
        sigma_data.append({"Substituent": "H",
        "sigma_p_minus": 0})
        return sigma_data

    def merge_with_energy_profiles(self, sigma_data):
        """
        Merges sigma_p_minus data with energy_profiles_data to get G_IM and G_TS2 for each substituent.
        Returns two lists for saving:
          - one with Substituent, sigma_p, G_IM
          - one with Substituent, sigma_p, G_TS2
        """
        energy_dict = {}
        ts2_dict = {}

        for entry in self.energy_profiles_data:
            sub = entry.get("Substituent")
            subdir = entry.get("Subdir")
            qhG_rel_kcal = entry.get("qhG rel (kcal/mol)")
            if sub is None or subdir is None or qhG_rel_kcal is None or qhG_rel_kcal == "":
                continue
            try:
                energy = float(qhG_rel_kcal)
            except ValueError:
                continue

            if sub not in energy_dict:
                energy_dict[sub] = {}
            energy_dict[sub][subdir] = energy

        for sub, stages in energy_dict.items():
            ts2 = stages.get("TS2")
            im = stages.get("IM")
            if ts2 is not None and im is not None:
                ts2_dict[sub] = ts2

        merged_G_IM = []
        merged_G_TS2 = []

        for item in sigma_data:
            sub = item["Substituent"]
            sigma_p_minus = item["sigma_p_minus"]
            g_im = energy_dict.get(sub, {}).get("IM")
            g_ts2 = ts2_dict.get(sub)

            merged_G_IM.append({
                "Substituent": sub,
                "sigma_p_minus": sigma_p_minus,
                "G_IM (kcal/mol)": g_im if g_im is not None else ""
            })

            merged_G_TS2.append({
                "Substituent": sub,
                "sigma_p_minus": sigma_p_minus,
                "G_TS2 (kcal/mol)": g_ts2 if g_ts2 is not None else ""
            })

        return merged_G_IM, merged_G_TS2

    def save_to_csv(self, merged_G_IM, merged_G_TS2):
        '''Creates Fig_3D.csv and Fig_3F.csv'''
        with open("Fig_3D.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["Substituent", "sigma_p_minus", "G_IM (kcal/mol)"])
            writer.writeheader()
            writer.writerows(merged_G_IM)
        print("\nCreated Fig_3D.csv\n")

        with open("Fig_3F.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["Substituent", "sigma_p_minus", "G_TS2 (kcal/mol)"])
            writer.writeheader()
            writer.writerows(merged_G_TS2)
        print("\nCreated Fig_3F.csv\n")

    def run(self):
        hammett_entries = self.parse_hammett_data()
        self.save_hammett_data_to_csv(hammett_entries)
        sigma_data = self.calculate_sigma_p(hammett_entries)
        merged_G_IM, merged_G_TS2 = self.merge_with_energy_profiles(sigma_data)
        self.save_to_csv(merged_G_IM, merged_G_TS2)


def extract_final_geometries():
    """
    Creates SI_structures_trj.xyz with all final Cartesian coordinates from directories Energy_profiles,
    E_HU, Hammett_constants by parsing .log files from those directories, in a specific order.
    """
    dirs = ["Energy_profiles", "Hammett_constants", "E_HU"]
    output_path = "SI_structures_trj.xyz"
    order_entry = ["3'", "TS1", "IM", "TS2", "4'", "simplified_IM"]

    print("\n---Creating _trj.xyz file with final structure Cartesian coordinates (angstroem)---\n")

    all_files = []
    for dir_name in dirs:
        base_dir = Utilities.get_dir(dir_name, "SI_archive")
        for root, _, files in os.walk(base_dir):
            for fname in files:
                if fname.endswith(".log"):
                    filepath = os.path.join(root, fname)

                    entry_type = None
                    for entry in order_entry:
                        if entry in fname or entry in os.path.basename(root):
                            entry_type = entry
                            break

                    substituent = None
                    for entry in order_entry:
                        if fname.startswith(entry + "_"):
                            substituent = fname[len(entry)+1:-4]
                            break
                    if substituent is None and "_" in fname:
                        substituent = fname.split("_", 1)[-1][:-4]
                    all_files.append({
                        "filepath": filepath,
                        "fname": fname,
                        "entry_type": entry_type,
                        "substituent": substituent
                    })

    ordered_files = []
    used_files = set()
    for entry in order_entry:

        entry_files = [f for f in all_files if f["entry_type"] == entry]
        entry_files.sort(key=lambda x: x["substituent"] or "")
        ordered_files.extend(entry_files)
        used_files.update(f["filepath"] for f in entry_files)

    remaining_files = [f for f in all_files if f["filepath"] not in used_files]
    remaining_files.sort(key=lambda x: (x["substituent"] or "", x["fname"]))
    ordered_files.extend(remaining_files)

    total_files = len(ordered_files)
    processed_files = 0

    with open(output_path, 'w', encoding='utf-8') as out_f:
        for fileinfo in ordered_files:
            filepath = fileinfo["filepath"]
            fname = fileinfo["fname"]
            try:
                with open(filepath, 'r', encoding='utf-8') as file:
                    lines = file.readlines()
                fe_indices = [i for i, line in enumerate(lines) if "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" in line]
                if not fe_indices:
                    print(f"  No FINAL ENERGY EVALUATION block found in {filepath}")
                    continue
                fe_idx = fe_indices[-1]
                cc_idx = None
                for i in range(fe_idx, min(fe_idx+100, len(lines))):
                    if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[i]:
                        cc_idx = i
                        break
                if cc_idx is None:
                    print(f"  No ANGSTROEM coordinate block found in {filepath}")
                    continue
                coord_start = cc_idx + 1
                while coord_start < len(lines) and lines[coord_start].strip().startswith('-'):
                    coord_start += 1
                coord_lines = []
                for line in lines[coord_start:]:
                    if line.strip().startswith('-') or "CARTESIAN COORDINATES" in line or "ATOMIC COORDINATES" in line or "END" in line:
                        break
                    if line.strip():
                        coord_lines.append(line.rstrip())
                if not coord_lines:
                    print(f"  No coordinates found in {filepath}")
                    continue
                base_name_no_ext = os.path.splitext(fname)[0]
                out_f.write(base_name_no_ext + "\n")
                for cl in coord_lines:
                    out_f.write(cl + "\n")
                out_f.write("\n")
                processed_files += 1
            except Exception as e:
                print(f"  Error processing {filepath}: {e}")

    print(f"\nTotal .log files found: {total_files}")
    print(f"Total files successfully processed: {processed_files}")
    print(f"Structures written to {output_path}\n")


def run_extraction():
    '''Runs extraction to create .csv files from computational data in ORCA .log files'''
    orbital_data_extractor = OrbitalEnergyExtractor()
    orbital_data = orbital_data_extractor.extract_orbital_data()
    energy_analyzer = EnergyProfileAnalyzer()
    energy_profile_data = energy_analyzer.run_analysis()
    orbital_data_extractor.save_to_csv(orbital_data, energy_profile_data) 
    hammett = HammettConstants(energy_profile_data)
    hammett.run()

def main():
    extract_final_geometries()
    run_extraction()

if __name__ == "__main__":
    main()

