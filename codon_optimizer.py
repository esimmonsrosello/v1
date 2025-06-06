#Need to fix the single optimisation option 6 chart output - its not outputting anything until you do batch option 6
#Need to sort out the format of option 2 in batch optimisation
#make the logo rounder like an app

import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, scrolledtext, font as tkFont
import threading
import platform
import logging
import re
import io
import base64
import time
import pandas as pd
from PIL import Image, ImageTk
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import tempfile

# --- Constants and Configuration ---
CODON_USAGE_TABLE = {
    'TTT': 0.58, 'TTC': 0.42, 'TTA': 0.14, 'TTG': 0.13, 'TAT': 0.59, 'TAC': 0.41, 'TAA': 0.61, 'TAG': 0.08, 'TGA': 0.31,
    'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'TGT': 0.46, 'TGC': 0.54, 'TGG': 1.00,
    'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47, 'CAT': 0.57, 'CAC': 0.43, 'CAA': 0.34, 'CAG': 0.66,
    'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49, 'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11,
    'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, 'ATG': 1.00, 'AAT': 0.49, 'AAC': 0.51, 'AAA': 0.74, 'AAG': 0.26,
    'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.25, 'AGT': 0.16, 'AGC': 0.25, 'AGA': 0.07, 'AGG': 0.04,
    'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35, 'GAT': 0.63, 'GAC': 0.37, 'GAA': 0.68, 'GAG': 0.32,
    'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, 'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15,
}

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AA_TO_CODONS = {}
for codon, aa in GENETIC_CODE.items():
    if aa not in AA_TO_CODONS: AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)

logger = None

class RichTextHandler(logging.Handler):
    def __init__(self, text_area):
        super().__init__()
        self.text_area = text_area
        self.text_area.tag_config("INFO", foreground="black")
        self.text_area.tag_config("WARNING", foreground="orange", font=(None, 9, "bold"))
        self.text_area.tag_config("ERROR", foreground="red", font=(None, 9, "bold"))
        self.text_area.tag_config("DEBUG", foreground="grey")

    def emit(self, record):
        msg = self.format(record)
        level_tag = record.levelname.upper()
        self.text_area.configure(state='normal')
        self.text_area.insert(tk.END, msg + '\n', (level_tag,))
        self.text_area.configure(state='disabled')
        self.text_area.see(tk.END)

def setup_logging(log_area_widget):
    global logger
    logger = logging.getLogger("CodonOptimizerApp")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    handler = RichTextHandler(log_area_widget)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False
    return logger

def validate_dna_sequence(sequence):
    return bool(re.fullmatch(r"[ATGC]+", sequence.upper()))

def get_optimized_codon(amino_acid, codon_table_prefs):
    possible_codons = AA_TO_CODONS.get(amino_acid, [])
    if not possible_codons: return None
    best_codon = None; max_freq = -1
    for codon in possible_codons:
        freq = codon_table_prefs.get(codon, 0)
        if freq > max_freq: max_freq = freq; best_codon = codon
    return best_codon if best_codon else (possible_codons[0] if possible_codons else None)

def optimize_dna_sequence(dna_sequence, codon_table):
    if not validate_dna_sequence(dna_sequence):
        raise ValueError("Invalid DNA: Non-ATGC characters found.")
    dna_sequence_upper = dna_sequence.upper()
    protein_sequence = ""
    optimized_dna = ""
    for i in range(0, len(dna_sequence_upper) - (len(dna_sequence_upper) % 3), 3):
        codon = dna_sequence_upper[i:i+3]
        protein_sequence += GENETIC_CODE.get(codon, "?")
    if protein_sequence.endswith("?"): protein_sequence = protein_sequence[:-1]
    for aa in protein_sequence:
        if aa == "*":
            stop_codons = [c for c in AA_TO_CODONS['*'] if c in codon_table]
            best_stop = max(stop_codons, key=lambda c: codon_table.get(c,0), default="TAA") if stop_codons else "TAA"
            optimized_dna += best_stop
            continue
        opt_codon = get_optimized_codon(aa, codon_table)
        if opt_codon: optimized_dna += opt_codon
        else:
            original_codon_for_aa = next((c for c in AA_TO_CODONS.get(aa,[])), None)
            if original_codon_for_aa: optimized_dna += original_codon_for_aa
            else: optimized_dna += "NNN"
    return optimized_dna, protein_sequence

def calculate_gc_content(dna_sequence):
    if not dna_sequence: return 0.0
    return (dna_sequence.upper().count('G') + dna_sequence.upper().count('C')) / len(dna_sequence) * 100

def count_codons(dna_sequence):
    counts = {}
    for i in range(0, len(dna_sequence) - (len(dna_sequence) % 3), 3):
        codon = dna_sequence[i:i+3].upper()
        counts[codon] = counts.get(codon, 0) + 1
    return counts

def analyze_plus1_frameshift_stops(dna_sequence):
    stop_counts = {'TAA': 0, 'TAG': 0, 'TGA': 0}
    shifted_sequence = dna_sequence[1:]
    for i in range(0, len(shifted_sequence) - (len(shifted_sequence) % 3), 3):
        codon = shifted_sequence[i:i+3].upper()
        if codon in stop_counts: stop_counts[codon] += 1
    return stop_counts

# --- Batch Processing Functions ---
def translate_dna(sequence):
    protein = ""
    for i in range(0, len(sequence) - (len(sequence) % 3), 3):
        codon = sequence[i:i+3].upper()
        protein += GENETIC_CODE.get(codon, "?")
    if protein.endswith("?") and len(sequence) % 3 != 0 :
        protein = protein[:-1]
    return protein

def get_codon_weights_row(sequence):
    codons_list = []
    weights = []
    for i in range(0, len(sequence) - (len(sequence) % 3), 3):
        codon = sequence[i:i+3].upper()
        codons_list.append(codon)
        weights.append(CODON_USAGE_TABLE.get(codon, 0.0))
    return weights, codons_list

def process_batch_sequences(logger_param, input_text, choice, filename_template, parent_window_for_messages):
    current_processing_logger = logger_param if logger_param else logger
    sequences_input = input_text.strip().split('>')
    results = []
    processed_names = set()

    for seq_block in sequences_input:
        if not seq_block.strip(): continue
        lines = seq_block.strip().split('\n')
        header = lines[0].strip()
        dna_sequence_parts = [line.strip().replace(" ", "") for line in lines[1:]]
        raw_dna_seq = "".join(dna_sequence_parts)
        name = header if header else f"UnnamedSequence_{len(processed_names)+1}"
        original_name = name; name_count = 1
        while name in processed_names: name = f"{original_name}_{name_count}"; name_count += 1
        processed_names.add(name)
        clean_seq = re.sub(r"[^ATGCatgc]", "", raw_dna_seq).upper()

        if not validate_dna_sequence(clean_seq):
            if current_processing_logger: current_processing_logger.warning(f"Seq {name} invalid. Skipping.")
            results.append({'Sequence_Name': name, 'Error': 'Invalid/empty DNA sequence.'})
            continue
        if len(clean_seq) < 3:
            if current_processing_logger: current_processing_logger.warning(f"Seq {name} too short. Skipping.")
            results.append({'Sequence_Name': name, 'Error': 'Sequence too short.'})
            continue
        try:
            protein_seq = translate_dna(clean_seq)
            result_data = {'Sequence_Name': name, 'Original_DNA': clean_seq, 'Protein': protein_seq}

            if choice == 1:
                optimized_dna, _ = optimize_dna_sequence(clean_seq, CODON_USAGE_TABLE)
                result_data.update({'Optimized_DNA': optimized_dna, 'Analysis_Type': 'Codon Optimization'})
                results.append(result_data)
            elif choice == 2: # Logic for choice 2 as per Turn 31 requirements
                weights, codons_list = get_codon_weights_row(clean_seq)
                amino_acids = [GENETIC_CODE.get(c, '?') for c in codons_list]
                for i in range(len(codons_list)):
                    codon_specific_data = {
                        'Sequence_Name': name,
                        'Analysis_Type': 'CAI Weight Analysis (Original)',
                        'Original_DNA': clean_seq,
                        'Protein': protein_seq,
                        'Position': i + 1,
                        'DNA_Codon': codons_list[i],
                        'CAI_Weight': weights[i],
                        'Amino_Acid': amino_acids[i]
                    }
                    results.append(codon_specific_data)
                continue
            elif choice == 3:
                gc_original = calculate_gc_content(clean_seq)
                result_data.update({'GC_Content_Original': f"{gc_original:.2f}%", 'Analysis_Type': 'GC Content (Original)'})
                results.append(result_data)
            elif choice == 4:
                codon_counts_dict = count_codons(clean_seq)
                result_data.update({'Original_Codon_Counts': str(codon_counts_dict), 'Analysis_Type': 'Codon Counts (Original)'})
                results.append(result_data)
            elif choice == 5:
                optimized_dna, _ = optimize_dna_sequence(clean_seq, CODON_USAGE_TABLE)
                gc_original = calculate_gc_content(clean_seq); gc_optimized = calculate_gc_content(optimized_dna)
                opt_weights, _ = get_codon_weights_row(optimized_dna)
                cai_opt_str = ','.join(f"{w:.4f}" for w in opt_weights) if opt_weights else ""
                result_data.update({'Optimized_DNA': optimized_dna, 'GC_Content_Original': f"{gc_original:.2f}%",
                                   'GC_Content_Optimized': f"{gc_optimized:.2f}%", 'CAI_Optimized': cai_opt_str,
                                   'Analysis_Type': 'Optimize & Analyze All'})
                results.append(result_data)
            elif choice == 6:
                optimized_dna, _ = optimize_dna_sequence(clean_seq, CODON_USAGE_TABLE)
                stops = analyze_plus1_frameshift_stops(optimized_dna)
                result_data.update({'Optimized_DNA': optimized_dna, 'Plus1_TAA_Count': stops.get('TAA',0),
                                   'Plus1_TAG_Count': stops.get('TAG',0), 'Plus1_TGA_Count': stops.get('TGA',0),
                                   'Plus1_Total_Stops': sum(stops.values()), 'Analysis_Type': '+1 Stops (Optimized)'})
                results.append(result_data)
        except Exception as e:
            if current_processing_logger: current_processing_logger.error(f"Error processing {name}: {e}", exc_info=True)
            results.append({'Sequence_Name': name, 'Error': str(e)})

    if not results:
        if current_processing_logger: current_processing_logger.info("No results generated.")
        messagebox.showinfo("Batch Info", "No results.", parent=parent_window_for_messages); return None

    all_keys_initial = set().union(*(d.keys() for d in results)); final_ordered_keys = []
    if choice == 6:
        preferred_order = ['Sequence_Name','Analysis_Type','Original_DNA','Protein','Optimized_DNA','Plus1_TAA_Count','Plus1_TAG_Count','Plus1_TGA_Count','Plus1_Total_Stops','Error']
        final_ordered_keys.extend([k for k in preferred_order if k in all_keys_initial])
        remaining_keys = sorted([k for k in all_keys_initial if k not in final_ordered_keys])
        final_ordered_keys.extend(remaining_keys)
    elif choice == 2: # Column ordering for choice 2 as per Turn 31 requirements
        preferred_order_choice2 = ['Sequence_Name','Analysis_Type','Original_DNA','Protein','Position','DNA_Codon','Amino_Acid','CAI_Weight','Error']
        actual_data_keys = set(results[0].keys()) if results and isinstance(results[0], dict) else all_keys_initial
        final_ordered_keys.extend([key for key in preferred_order_choice2 if key in actual_data_keys])
        # Add Error key if present and not already included, then other remaining keys
        if 'Error' in actual_data_keys and 'Error' not in final_ordered_keys: final_ordered_keys.append('Error')
        final_ordered_keys.extend(sorted([key for key in actual_data_keys if key not in final_ordered_keys]))
    else:
        final_ordered_keys = sorted(list(all_keys_initial)) # Original default for other choices

    standardized_results = [{k: row.get(k) for k in final_ordered_keys} for row in results]
    try:
        df = pd.DataFrame(standardized_results)
        if choice == 2 and 'CAI_Weight' in df.columns:
             df['CAI_Weight'] = df['CAI_Weight'].apply(lambda x: f"{x:.4f}" if isinstance(x, float) else x)
        ts = time.strftime("%Y%m%d-%H%M%S")
        out_fn = f"{filename_template}_batch_op{choice}_{ts}.xlsx"
        df.to_excel(out_fn, index=False, engine='xlsxwriter')
        if current_processing_logger: current_processing_logger.info(f"Batch complete. Saved to {out_fn}")
        messagebox.showinfo("Batch Complete", f"Finished.\nSaved to: {out_fn}", parent=parent_window_for_messages)
        return out_fn
    except Exception as e:
        if current_processing_logger: current_processing_logger.error(f"Batch save to Excel failed: {e}", exc_info=True)
        messagebox.showerror("Batch Save Error", f"Failed to save: {e}", parent=parent_window_for_messages)
        return None
# --- End of Batch Processing ---

# --- File Handling & Output --- (Original from Turn 39 script for single sequence ops)
def save_output_to_excel(data_dict_for_excel, filename_template, chart_image_path=None):
    global logger
    try:
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        filename = f"{filename_template}_{timestamp}.xlsx"
        df_data_sheets = {}
        if 'summary_stats' in data_dict_for_excel:
            if isinstance(data_dict_for_excel['summary_stats'], pd.DataFrame): df_data_sheets['Summary Statistics'] = data_dict_for_excel['summary_stats']
            elif isinstance(data_dict_for_excel['summary_stats'], dict): df_data_sheets['Summary Statistics'] = pd.DataFrame([data_dict_for_excel['summary_stats']])
        if 'sequences' in data_dict_for_excel:
            if isinstance(data_dict_for_excel['sequences'], pd.DataFrame): df_data_sheets['Sequences'] = data_dict_for_excel['sequences']
            elif isinstance(data_dict_for_excel['sequences'], dict): df_data_sheets['Sequences'] = pd.DataFrame([data_dict_for_excel['sequences']])
        if 'original_codon_counts' in data_dict_for_excel and data_dict_for_excel['original_codon_counts']:
            s = pd.Series(data_dict_for_excel['original_codon_counts'], name="Count"); s.index.name = "Codon"; df_data_sheets['Original Codon Counts'] = s.reset_index()
        if 'optimized_codon_counts' in data_dict_for_excel and data_dict_for_excel['optimized_codon_counts']:
            s = pd.Series(data_dict_for_excel['optimized_codon_counts'], name="Count"); s.index.name = "Codon"; df_data_sheets['Optimized Codon Counts'] = s.reset_index()
        if 'plus1_stop_counts' in data_dict_for_excel and data_dict_for_excel['plus1_stop_counts']:
            s = pd.Series(data_dict_for_excel['plus1_stop_counts'], name="Count"); s.index.name = "Stop_Codon_Type"; df_data_sheets['+1 Frame Stops'] = s.reset_index()
        with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
            for sheet_name, df in df_data_sheets.items():
                if not df.empty: df.to_excel(writer, sheet_name=sheet_name, index=False)
            if chart_image_path and os.path.exists(chart_image_path):
                if '+1 Frame Stops' in writer.sheets: writer.sheets['+1 Frame Stops'].insert_image('E2', chart_image_path)
                elif logger: logger.warning("Sheet '+1 Frame Stops' not found for chart.")
        if logger: logger.info(f"Output saved to {filename}")
        return filename
    except Exception as e:
        if logger: logger.error(f"Excel save failed: {e}", exc_info=True)
        messagebox.showerror("Save Error", f"Failed to save Excel: {e}", parent=None)
        return None

MACOS_ICON_B64 = """
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
"""

def set_macos_app_icon(root_window, logger_instance_param):
    current_logger = logger_instance_param if logger_instance_param else logger
    if platform.system() == "Darwin":
        try:
            if not hasattr(set_macos_app_icon, 'icon_image'):
                icon_data = base64.b64decode(MACOS_ICON_B64)
                set_macos_app_icon.icon_image = tk.PhotoImage(data=icon_data)
            if set_macos_app_icon.icon_image:
                root_window.iconphoto(True, set_macos_app_icon.icon_image)
                if current_logger: current_logger.info("App icon set for macOS.")
            elif current_logger: current_logger.warning("Icon image not loaded.")
        except Exception as e:
            if current_logger: current_logger.error(f"Setting macOS icon failed: {e}")

def create_gui(root):
    global logger
    root.title("Codon Optimizer Pro")
    root.geometry("900x700")
    main_paned_window = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
    main_paned_window.pack(fill=tk.BOTH, expand=True)
    left_frame_outer = ttk.Frame(main_paned_window, padding=10)
    main_paned_window.add(left_frame_outer, weight=35)
    left_frame = ttk.Frame(left_frame_outer)
    left_frame.pack(fill=tk.BOTH, expand=True)
    ttk.Label(left_frame, text="Enter DNA Sequence(s) (FASTA or raw):", font=(None, 11, "bold")).pack(pady=(0,5), anchor='w')
    input_sequence_text = scrolledtext.ScrolledText(left_frame, height=12, width=50, relief=tk.SOLID, borderwidth=1)
    input_sequence_text.pack(fill=tk.BOTH, expand=True)
    input_sequence_text.focus()
    ttk.Label(left_frame, text="Select Operation:", font=(None, 11, "bold")).pack(pady=(10,5), anchor='w')
    operation_var = tk.IntVar(value=1)
    choices = [("Optimize Sequence",1), ("Calculate GC Content",2), ("Count Codons (Original)",3),
               ("Translate to Protein",4), ("Optimize & Analyze All",5), ("Analyze +1 Frameshift Stops (Optimized)",6)]
    for txt, val in choices: ttk.Radiobutton(left_frame, text=txt, variable=operation_var, value=val).pack(anchor='w')
    combine_outputs_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(left_frame, text="Accumulate results from multiple runs", variable=combine_outputs_var).pack(pady=(5,0), anchor='w')
    ttk.Label(left_frame, text="Filename Template for Save:", font=(None, 10, "bold")).pack(pady=(10,2), anchor='w')
    filename_template_var = tk.StringVar(value="optimization_output")
    ttk.Entry(left_frame, textvariable=filename_template_var, width=35).pack(anchor='w', fill=tk.X, expand=True)
    right_frame_outer = ttk.Frame(main_paned_window, padding=10)
    main_paned_window.add(right_frame_outer, weight=65)
    right_frame_notebook = ttk.Notebook(right_frame_outer)
    right_frame_notebook.pack(fill=tk.BOTH, expand=True)
    output_tab = ttk.Frame(right_frame_notebook); right_frame_notebook.add(output_tab, text='Results / Analysis')
    output_sequence_text = scrolledtext.ScrolledText(output_tab, height=15, width=60, relief=tk.SOLID, borderwidth=1, state='disabled', wrap=tk.WORD)
    output_sequence_text.pack(fill=tk.BOTH, expand=True, pady=(0,5))
    chart_frame = ttk.Frame(output_tab)
    chart_frame.pack(fill=tk.BOTH, expand=True, pady=(5,0))
    ttk.Label(chart_frame, text="Chart (if any) will be displayed here or in a new window.").pack(expand=True)
    log_tab = ttk.Frame(right_frame_notebook); right_frame_notebook.add(log_tab, text='Log')
    log_area_text = scrolledtext.ScrolledText(log_tab, height=6, width=60, relief=tk.SOLID, borderwidth=1, state='disabled', wrap=tk.WORD)
    log_area_text.pack(fill=tk.BOTH, expand=True)
    logger = setup_logging(log_area_text)
    button_frame = ttk.Frame(left_frame)
    button_frame.pack(fill=tk.X, pady=(15,0), side=tk.BOTTOM)
    run_button = ttk.Button(button_frame, text="Run", command=lambda: run_optimization_thread(
        root, logger, input_sequence_text, output_sequence_text, operation_var,
        filename_template_var, combine_outputs_var, chart_frame))
    run_button.pack(side=tk.LEFT, padx=(0,10), fill=tk.X, expand=True)
    ttk.Button(button_frame, text="Clear Inputs", command=lambda: input_sequence_text.delete(1.0, tk.END)).pack(side=tk.LEFT, fill=tk.X, expand=True)
    if not hasattr(root, 'combined_data_storage'):
        root.combined_data_storage = {"summary_stats": pd.DataFrame(), "original_codon_counts": {}, "optimized_codon_counts": {}, "plus1_stop_counts": {}, "sequences": pd.DataFrame()}
    set_macos_app_icon(root, logger)

def run_optimization_thread(root, logger_param, input_widget, output_widget, operation_var,
                            filename_template_var, combine_outputs_var, chart_display_frame):
    run_button = None
    try:
        button_frame_parent = input_widget.master
        for child in button_frame_parent.winfo_children():
            if isinstance(child, ttk.Frame):
                for btn_candidate in child.winfo_children():
                    if isinstance(btn_candidate, ttk.Button) and "Run" in btn_candidate.cget("text"):
                        run_button = btn_candidate; break
                if run_button: break
    except Exception:
        if logger_param: logger_param.debug("Run button not found to disable.")
    if run_button: run_button.config(state=tk.DISABLED)
    thread = threading.Thread(target=run_optimization, args=(
        root, logger_param, input_widget, output_widget, operation_var,
        filename_template_var, combine_outputs_var, chart_display_frame, run_button))
    thread.daemon = True; thread.start()

current_temp_chart_file = None

def run_optimization(root, logger_instance, input_widget, output_widget, operation_var,
                     filename_template_var, combine_outputs_var, chart_display_frame, run_button_widget):
    global current_temp_chart_file
    effective_logger = logger_instance if logger_instance else logger
    parent_window = root
    try:
        input_text = input_widget.get(1.0, tk.END).strip()
        choice = operation_var.get()
        filename_template = filename_template_var.get().strip() or "optimization_output"
        actual_dna_sequence = input_text
        if input_text.startswith(">"):
            try:
                header, *seq_lines = input_text.split('\n')
                actual_dna_sequence = "".join(s.strip() for s in seq_lines)
                if effective_logger: effective_logger.info(f"Processing FASTA: {header.strip()}")
            except Exception as e_fasta:
                if effective_logger: effective_logger.error(f"FASTA parsing error: {e_fasta}")
                messagebox.showerror("Input Error", "Error parsing FASTA.", parent=parent_window); return
        actual_dna_sequence = re.sub(r"[^ATGCatgc]", "", actual_dna_sequence).upper()
        is_combined_op_on_empty = choice == 5 and combine_outputs_var.get() and root.combined_data_storage['summary_stats'].empty() and not actual_dna_sequence
        if not actual_dna_sequence and not is_combined_op_on_empty :
            messagebox.showwarning("Input Error", "Please enter DNA sequence.", parent=parent_window)
            if effective_logger: effective_logger.warning("No DNA input."); return
        if not validate_dna_sequence(actual_dna_sequence) and not is_combined_op_on_empty and choice != 5 :
             messagebox.showerror("Input Error", "Invalid DNA. Use A,T,G,C only.", parent=parent_window)
             if effective_logger: effective_logger.error("Invalid DNA sequence."); return
        output_widget.configure(state='normal'); output_widget.delete(1.0, tk.END); output_widget.configure(state='disabled')
        for widget in chart_display_frame.winfo_children(): widget.destroy()
        ttk.Label(chart_display_frame, text="Chart (if any) will appear here or new window.").pack(expand=True)
        if effective_logger: effective_logger.info(f"Op: {choice}, Accumulate: {combine_outputs_var.get()}")
        current_run_data = {"summary_stats": pd.DataFrame(), "original_codon_counts": {}, "optimized_codon_counts": {}, "plus1_stop_counts": {}, "sequences": pd.DataFrame()}
        chart_image_to_embed_in_excel = None
        if current_temp_chart_file and os.path.exists(current_temp_chart_file):
            try: os.remove(current_temp_chart_file); effective_logger.info(f"Removed old temp chart: {current_temp_chart_file}")
            except Exception: pass
            current_temp_chart_file = None
        plus1_stop_counts_data = {}
        if choice == 1 or choice == 5:
            optimized_dna, protein_seq = optimize_dna_sequence(actual_dna_sequence, CODON_USAGE_TABLE)
            current_run_data["sequences"] = pd.DataFrame({'Original Sequence': [actual_dna_sequence], 'Optimized Sequence': [optimized_dna], 'Protein Sequence': [protein_seq]})
            current_run_data["summary_stats"] = pd.DataFrame({'Metric': ['Original GC%', 'Optimized GC%', 'Original Length', 'Optimized Length'], 'Value': [f"{calculate_gc_content(actual_dna_sequence):.2f}%", f"{calculate_gc_content(optimized_dna):.2f}%", len(actual_dna_sequence), len(optimized_dna)]})
            current_run_data["original_codon_counts"] = count_codons(actual_dna_sequence)
            current_run_data["optimized_codon_counts"] = count_codons(optimized_dna)
            if choice == 1: _update_output_text(output_widget, f"Original:\n{actual_dna_sequence}\n\nOptimized:\n{optimized_dna}\n\nProtein:\n{protein_seq}\n\n" + current_run_data["summary_stats"].to_string(index=False))
        if choice == 2 or choice == 5:
            seq_for_gc = actual_dna_sequence
            seq_type_label = "Original"
            if choice == 5 and 'Optimized Sequence' in current_run_data["sequences"].columns and not current_run_data["sequences"]['Optimized Sequence'].empty:
                seq_for_gc = current_run_data["sequences"]['Optimized Sequence'].iloc[0]; seq_type_label = "Optimized"
            gc_val = calculate_gc_content(seq_for_gc)
            gc_df_row = pd.DataFrame({'Metric': [f'{seq_type_label} GC%'], 'Value': [f"{gc_val:.2f}%"]})
            if choice == 2: _update_output_text(output_widget, f"{seq_type_label} GC Content: {gc_val:.2f}%"); current_run_data["summary_stats"] = gc_df_row
            elif not current_run_data["summary_stats"].empty: current_run_data["summary_stats"] = pd.concat([current_run_data["summary_stats"], gc_df_row], ignore_index=True)
            else: current_run_data["summary_stats"] = gc_df_row
        if choice == 3 or choice == 5:
            counts = count_codons(actual_dna_sequence)
            if choice == 3: _update_output_text(output_widget, "Original Codon Counts:\n" + "\n".join([f"{c}: {v}" for c,v in counts.items()]))
            current_run_data["original_codon_counts"] = counts
        if choice == 4 or choice == 5:
            protein_s = optimize_dna_sequence(actual_dna_sequence, CODON_USAGE_TABLE)[1] if choice == 4 else current_run_data["sequences"]['Protein Sequence'].iloc[0]
            if choice == 4: _update_output_text(output_widget, f"Protein Sequence:\n{protein_s}")
        if choice == 6 or choice == 5:
            target_seq_for_plus1 = actual_dna_sequence
            if choice == 5 and 'Optimized Sequence' in current_run_data["sequences"].columns and not current_run_data["sequences"]['Optimized Sequence'].empty:
                 target_seq_for_plus1 = current_run_data["sequences"]['Optimized Sequence'].iloc[0]
            elif choice == 6:
                 target_seq_for_plus1, _ = optimize_dna_sequence(actual_dna_sequence, CODON_USAGE_TABLE)
            plus1_stop_counts_data = analyze_plus1_frameshift_stops(target_seq_for_plus1)
            current_run_data["plus1_stop_counts"] = plus1_stop_counts_data
            if choice == 6: _update_output_text(output_widget, "+1 Stops (Optimized):\n" + "\n".join([f"{c}: {v}" for c,v in plus1_stop_counts_data.items()]))
        if choice == 5:
            summary_text = current_run_data["summary_stats"].to_string(index=False) if not current_run_data["summary_stats"].empty else "N/A"
            disp5 = (f"--- Optimized Sequence & Analysis ---\n"
                     f"Original:\n{current_run_data['sequences']['Original Sequence'].iloc[0]}\n\n"
                     f"Optimized:\n{current_run_data['sequences']['Optimized Sequence'].iloc[0]}\n\n"
                     f"Protein:\n{current_run_data['sequences']['Protein Sequence'].iloc[0]}\n\n"
                     f"--- Summary Stats ---\n{summary_text}\n\n"
                     f"--- Original Codon Counts ---\n{', '.join(f'{k}:{v}' for k,v in current_run_data['original_codon_counts'].items())}\n\n"
                     f"--- Optimized Codon Counts ---\n{', '.join(f'{k}:{v}' for k,v in current_run_data['optimized_codon_counts'].items())}\n\n"
                     f"--- +1 Stops (Optimized) ---\n{', '.join(f'{k}:{v}' for k,v in current_run_data['plus1_stop_counts'].items())}")
            _update_output_text(output_widget, disp5)
            counts_chart5 = current_run_data["plus1_stop_counts"]
            labels_all5 = ['TAA','TAG','TGA']; sizes_all5 = [counts_chart5.get(lbl,0) for lbl in labels_all5]
            valid_lbl5 = [l for i,l in enumerate(labels_all5) if sizes_all5[i]>0]; valid_sz5 = [s for s in sizes_all5 if s>0]
            if sum(valid_sz5) > 0:
                try:
                    tfd5, tpath5 = tempfile.mkstemp(suffix=".png",prefix="chart_"); os.close(tfd5)
                    current_temp_chart_file = tpath5; chart_image_to_embed_in_excel = current_temp_chart_file
                    fig5,ax5=plt.subplots(); wed5,_,aut5=ax5.pie(valid_sz5,colors=[{'TAA':'#111d94','TAG':'#752db0','TGA':'#b34c9d'}[lbl] for lbl in valid_lbl5],autopct='%1.1f%%',startangle=90,textprops={'fontsize':9})
                    ax5.legend(wed5,valid_lbl5,title="+1 Stops(Opt)",loc="center left",bbox_to_anchor=(1,0.5),fontsize=8); ax5.axis('equal'); plt.tight_layout()
                    plt.savefig(chart_image_to_embed_in_excel,bbox_inches='tight'); plt.close(fig5)
                    if effective_logger: effective_logger.info(f"Chart for choice 5 to {chart_image_to_embed_in_excel}")
                    img5=Image.open(chart_image_to_embed_in_excel); img5.thumbnail((chart_frame.winfo_width() or 300, chart_frame.winfo_height() or 200),Image.Resampling.LANCZOS)
                    pho5=ImageTk.PhotoImage(img5); lbl_wid5=ttk.Label(chart_frame,image=pho5); lbl_wid5.image=pho5; lbl_wid5.pack(expand=True)
                except Exception as e_c5: effective_logger.error(f"Chart gen/disp choice 5 fail:{e_c5}"); chart_image_to_embed_in_excel=None; current_temp_chart_file=None
            else: effective_logger.info("No +1 stops for choice 5 chart."); chart_image_to_embed_in_excel=None; current_temp_chart_file=None
        data_to_proc = {k:v for k,v in current_run_data.items() if (isinstance(v,pd.DataFrame) and not v.empty) or (isinstance(v,dict) and v)}
        if combine_outputs_var.get():
            if "summary_stats" in data_to_proc: root.combined_data_storage["summary_stats"]=pd.concat([root.combined_data_storage["summary_stats"],data_to_proc["summary_stats"]]).drop_duplicates().reset_index(drop=True)
            for k in ["original_codon_counts","optimized_codon_counts","plus1_stop_counts"]:
                if k in data_to_proc:
                    for cod,cnt in data_to_proc[k].items(): root.combined_data_storage[k][cod]=root.combined_data_storage[k].get(cod,0)+cnt
            if "sequences" in data_to_proc: root.combined_data_storage["sequences"]=pd.concat([root.combined_data_storage["sequences"],data_to_proc["sequences"]]).reset_index(drop=True)
            if effective_logger: effective_logger.info("Data accumulated.")
            if choice in [1,5,6] and any((isinstance(v,pd.DataFrame)and not v.empty)or(isinstance(v,dict)and v) for v in root.combined_data_storage.values()):
                out_fn_cmb=save_output_to_excel(root.combined_data_storage,f"{filename_template}_COMBINED",chart_image_path=current_temp_chart_file)
                if out_fn_cmb: messagebox.showinfo("Complete (Combined)",f"Accumulated data to: {out_fn_cmb}",parent=parent_window)
            if root and root.winfo_exists(): set_macos_app_icon(root,effective_logger)
        else: # Not accumulating - This block contains the Turn 17/18 fixes for choice 6
            out_fn_indiv=None; excel_chart_path=chart_image_to_embed_in_excel
            if choice == 6:
                lbls6=['TAA','TAG','TGA']; siz6=[plus1_stop_counts_data.get(lbl,0) for lbl in lbls6]
                v_lbl6=[l for i,l in enumerate(lbls6) if siz6[i]>0]; v_siz6=[s for s in siz6 if s>0]
                cmap6={'TAA':'#111d94','TAG':'#752db0','TGA':'#b34c9d'}; v_col6=[cmap6[lbl] for lbl in v_lbl6 if lbl in cmap6]
                if sum(v_siz6)>0:
                    try:
                        fig6,ax6=plt.subplots()
                        wedges6, _, autotexts6 = ax6.pie(v_siz6,colors=v_col6,autopct='%1.1f%%',startangle=90,textprops={'fontsize':9}) # Has autopct
                        ax6.legend(wedges6,v_lbl6,title="+1 Stop Codons",loc="center left",bbox_to_anchor=(1,0.5),fontsize=8,title_fontsize=9) # Has v_lbl6 (valid_chart_labels)
                        ax6.axis('equal'); plt.tight_layout(); plt.show(block=False) # Has block=False
                        if effective_logger: effective_logger.info("Chart choice 6 (single) non-blocking.")
                    except Exception as e6: effective_logger.error(f"Chart display choice 6 failed:{e6}")
                else: effective_logger.info("No +1 stops for choice 6 chart (single).")
                excel_data6={'plus1_stop_counts':plus1_stop_counts_data}
                if actual_dna_sequence: excel_data6['sequences']=pd.DataFrame({'Original Sequence':[actual_dna_sequence],'Protein Sequence':[optimize_dna_sequence(actual_dna_sequence,CODON_USAGE_TABLE)[1]]})
                out_fn_indiv=save_output_to_excel(excel_data6,filename_template,chart_image_path=None); excel_chart_path=None
            elif data_to_proc: out_fn_indiv=save_output_to_excel(data_to_proc,filename_template,chart_image_path=excel_chart_path)
            else: effective_logger.info("No data to save.")
            if out_fn_indiv: messagebox.showinfo("Complete",f"Op completed.\nSaved to: {out_fn_indiv}",parent=parent_window)
            if root and root.winfo_exists(): set_macos_app_icon(root,effective_logger) # This is the corrected location
            if excel_chart_path and os.path.exists(excel_chart_path):
                try: os.remove(excel_chart_path); effective_logger.info(f"Temp chart {excel_chart_path} removed."); current_temp_chart_file=None
                except Exception as e_rem: effective_logger.error(f"Failed to remove temp chart {excel_chart_path}:{e_rem}")
    except ValueError as ve:
        if effective_logger: effective_logger.error(f"ValueError:{ve}",exc_info=True)
        messagebox.showerror("Processing Error",str(ve),parent=parent_window)
    except Exception as e:
        if effective_logger: effective_logger.critical(f"Critical error:{e}",exc_info=True)
        messagebox.showerror("Unexpected Error",f"Error:{e}",parent=parent_window)
    finally:
        if run_button_widget and isinstance(run_button_widget,ttk.Button):
            try:
                if run_button_widget.winfo_exists(): run_button_widget.config(state=tk.NORMAL)
            except tk.TclError:
                if effective_logger: effective_logger.warning("Run button widget no longer exists.")

def _update_output_text(output_widget, text):
    output_widget.configure(state='normal'); output_widget.delete(1.0,tk.END)
    output_widget.insert(tk.END,text); output_widget.configure(state='disabled'); output_widget.see(tk.END)

def main():
    global logger; root=tk.Tk()
    style=ttk.Style(); style.theme_use('clam')
    style.configure("TButton",padding=5,relief="flat",font=(None,10)); style.map("TButton",background=[('active','#c0c0c0')])
    logging.basicConfig(level=logging.INFO); logger=logging.getLogger("CodonOptimizerApp_Main")
    create_gui(root)
    default_font=tkFont.nametofont("TkDefaultFont"); default_font.configure(size=10)
    root.option_add("*Font",default_font)
    text_font=tkFont.nametofont("TkTextFont"); text_font.configure(size=10)
    def on_closing():
        if logger: logger.info("App closing."); root.destroy()
    root.protocol("WM_DELETE_WINDOW",on_closing); root.mainloop()

if __name__ == "__main__":
    main()
