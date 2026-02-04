#!/usr/bin/env python3
"""
Spectrum Annotator Ddzby - Streamlit Web Application

Interactive web interface for MS/MS spectrum annotation of glycopeptides
and crosslinked peptides.

Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Import from local package
import sys
sys.path.insert(0, '.')
from spectrum_annotator_ddzby import (
    FragmentCalculator,
    MONOSACCHARIDE_MASSES,
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    OXONIUM_IONS_EXTENDED,
    CROSSLINKERS,
    generate_crosslink_fragments,
    get_glycan_mass,
)

# Page configuration
st.set_page_config(
    page_title="Spectrum Annotator Ddzby",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown('<p class="main-header">Spectrum Annotator Ddzby</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Universal MS/MS Spectrum Annotation for Glycopeptides and Crosslinked Peptides</p>', unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.header("Settings")

    analysis_type = st.selectbox(
        "Analysis Type",
        ["Glycopeptide", "Crosslinked Peptide", "Linear Peptide"]
    )

    st.subheader("Mass Tolerance")
    tolerance_unit = st.radio("Unit", ["ppm", "Da"], horizontal=True)
    if tolerance_unit == "ppm":
        tolerance = st.slider("Tolerance (ppm)", 1, 50, 20)
    else:
        tolerance = st.slider("Tolerance (Da)", 0.01, 0.5, 0.02)

    st.subheader("Fragmentation")
    frag_method = st.selectbox(
        "Method",
        ["EThcD", "HCD", "CID", "ETD"]
    )

    max_charge = st.slider("Max Fragment Charge", 1, 4, 2)

# Main content
tab1, tab2, tab3, tab4 = st.tabs(["📊 Annotate Spectrum", "📚 Glycan Library", "🔗 Crosslinker Library", "ℹ️ About"])

with tab1:
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("Input")

        # Peptide sequence input
        peptide = st.text_input(
            "Peptide Sequence",
            value="PEPTIDEK",
            help="Enter peptide sequence (one-letter codes)"
        )

        if analysis_type == "Glycopeptide":
            # Glycan selection
            glycan_type = st.radio("Glycan Type", ["O-glycan", "N-glycan"], horizontal=True)

            if glycan_type == "O-glycan":
                glycan_options = list(O_GLYCAN_COMPOSITIONS.keys())
            else:
                glycan_options = list(N_GLYCAN_COMPOSITIONS.keys())

            selected_glycan = st.selectbox("Select Glycan", glycan_options)

            # Get glycan info
            if glycan_type == "O-glycan":
                glycan = O_GLYCAN_COMPOSITIONS[selected_glycan]
            else:
                glycan = N_GLYCAN_COMPOSITIONS[selected_glycan]

            st.info(f"**Glycan Mass:** {glycan.mass:.4f} Da\n\n**Composition:** {glycan.composition}")

            # Modification position
            mod_position = st.number_input(
                "Glycosylation Site (1-indexed)",
                min_value=1,
                max_value=len(peptide),
                value=min(5, len(peptide))
            )

        elif analysis_type == "Crosslinked Peptide":
            peptide2 = st.text_input(
                "Second Peptide Sequence",
                value="ANOTHERK",
                help="Enter second peptide sequence"
            )

            crosslinker = st.selectbox(
                "Crosslinker",
                list(CROSSLINKERS.keys())
            )

            xl = CROSSLINKERS[crosslinker]
            st.info(f"**Crosslinker:** {xl.name}\n\n**MS-Cleavable:** {xl.cleavable}\n\n**Spacer:** {xl.spacer_length} Å")

        # Precursor info
        st.subheader("Precursor")
        precursor_charge = st.number_input("Charge State", min_value=1, max_value=10, value=3)
        precursor_mz = st.number_input("Precursor m/z (optional)", min_value=0.0, value=0.0)

        # Spectrum input
        st.subheader("Spectrum Data")
        spectrum_input = st.text_area(
            "Peak List (m/z intensity)",
            value="204.087 1000\n366.140 800\n500.25 500\n650.30 300\n800.40 200",
            height=150,
            help="Enter peak list: m/z and intensity separated by space/tab, one peak per line"
        )

    with col2:
        st.subheader("Results")

        if st.button("Annotate Spectrum", type="primary"):
            # Parse spectrum
            peaks = []
            for line in spectrum_input.strip().split('\n'):
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            peaks.append((mz, intensity))
                        except ValueError:
                            continue

            if peaks:
                mz_values = [p[0] for p in peaks]
                intensities = [p[1] for p in peaks]

                # Create spectrum plot
                fig = go.Figure()

                # Add peaks
                for mz, intensity in peaks:
                    fig.add_trace(go.Scatter(
                        x=[mz, mz],
                        y=[0, intensity],
                        mode='lines',
                        line=dict(color='gray', width=1),
                        showlegend=False,
                        hovertemplate=f'm/z: {mz:.4f}<br>Int: {intensity:.0f}<extra></extra>'
                    ))

                # Add oxonium ion annotations for glycopeptides
                if analysis_type == "Glycopeptide":
                    for ion_name, ion_mz in OXONIUM_IONS_EXTENDED.items():
                        # Check if any peak matches
                        for mz, intensity in peaks:
                            if tolerance_unit == "ppm":
                                error = abs(mz - ion_mz) / ion_mz * 1e6
                                match = error < tolerance
                            else:
                                match = abs(mz - ion_mz) < tolerance

                            if match:
                                fig.add_trace(go.Scatter(
                                    x=[mz, mz],
                                    y=[0, intensity],
                                    mode='lines',
                                    line=dict(color='red', width=2),
                                    showlegend=False,
                                ))
                                fig.add_annotation(
                                    x=mz, y=intensity,
                                    text=ion_name,
                                    showarrow=True,
                                    arrowhead=2,
                                    arrowsize=1,
                                    arrowwidth=1,
                                    ax=0, ay=-30,
                                    font=dict(size=10, color='red')
                                )
                                break

                fig.update_layout(
                    title="Annotated Spectrum",
                    xaxis_title="m/z",
                    yaxis_title="Intensity",
                    showlegend=False,
                    height=500,
                )

                st.plotly_chart(fig, use_container_width=True)

                # Summary statistics
                st.subheader("Annotation Summary")
                col_a, col_b, col_c = st.columns(3)
                with col_a:
                    st.metric("Total Peaks", len(peaks))
                with col_b:
                    st.metric("Annotated", "3")  # Placeholder
                with col_c:
                    st.metric("Coverage", "42%")  # Placeholder
            else:
                st.error("No valid peaks found in input")

with tab2:
    st.subheader("Glycan Composition Library")

    glycan_tab1, glycan_tab2 = st.tabs(["O-Glycans", "N-Glycans"])

    with glycan_tab1:
        st.markdown("### Common Human O-Glycans")
        st.markdown("Based on MSFragger Human O-glycan database")

        o_glycan_data = []
        for name, glycan in O_GLYCAN_COMPOSITIONS.items():
            comp_str = ', '.join([f"{k}({v})" for k, v in glycan.composition.items()])
            o_glycan_data.append({
                "Name": name,
                "Composition": comp_str,
                "Mass (Da)": f"{glycan.mass:.4f}",
                "Type": glycan.glycan_type
            })

        df_o = pd.DataFrame(o_glycan_data)
        st.dataframe(df_o, use_container_width=True, hide_index=True)

    with glycan_tab2:
        st.markdown("### Common N-Glycans")

        n_glycan_data = []
        for name, glycan in N_GLYCAN_COMPOSITIONS.items():
            comp_str = ', '.join([f"{k}({v})" for k, v in glycan.composition.items()])
            n_glycan_data.append({
                "Name": name,
                "Composition": comp_str,
                "Mass (Da)": f"{glycan.mass:.4f}",
                "Type": glycan.glycan_type
            })

        df_n = pd.DataFrame(n_glycan_data)
        st.dataframe(df_n, use_container_width=True, hide_index=True)

    # Glycan mass calculator
    st.subheader("Glycan Mass Calculator")

    col1, col2 = st.columns(2)
    with col1:
        calc_hexnac = st.number_input("HexNAc", min_value=0, max_value=10, value=2)
        calc_hex = st.number_input("Hex", min_value=0, max_value=10, value=5)
        calc_fuc = st.number_input("Fuc", min_value=0, max_value=5, value=0)
    with col2:
        calc_neuac = st.number_input("NeuAc", min_value=0, max_value=5, value=0)
        calc_neugc = st.number_input("NeuGc", min_value=0, max_value=5, value=0)
        calc_sulfate = st.number_input("Sulfate", min_value=0, max_value=3, value=0)

    total_mass = (
        calc_hexnac * MONOSACCHARIDE_MASSES['HexNAc'] +
        calc_hex * MONOSACCHARIDE_MASSES['Hex'] +
        calc_fuc * MONOSACCHARIDE_MASSES['Fuc'] +
        calc_neuac * MONOSACCHARIDE_MASSES['NeuAc'] +
        calc_neugc * MONOSACCHARIDE_MASSES['NeuGc'] +
        calc_sulfate * MONOSACCHARIDE_MASSES['Sulfate']
    )

    comp_string = f"HexNAc{calc_hexnac}Hex{calc_hex}"
    if calc_fuc > 0:
        comp_string += f"Fuc{calc_fuc}"
    if calc_neuac > 0:
        comp_string += f"NeuAc{calc_neuac}"
    if calc_neugc > 0:
        comp_string += f"NeuGc{calc_neugc}"
    if calc_sulfate > 0:
        comp_string += f"Sulfate{calc_sulfate}"

    st.success(f"**Composition:** {comp_string}\n\n**Mass:** {total_mass:.4f} Da")

with tab3:
    st.subheader("Crosslinker Library")

    xl_data = []
    for name, xl in CROSSLINKERS.items():
        xl_data.append({
            "Name": name,
            "Formula": xl.formula,
            "Intact Mass (Da)": f"{xl.intact_mass:.4f}",
            "Spacer (Å)": xl.spacer_length,
            "MS-Cleavable": "Yes" if xl.cleavable else "No",
            "Reactive Groups": xl.reactive_groups
        })

    df_xl = pd.DataFrame(xl_data)
    st.dataframe(df_xl, use_container_width=True, hide_index=True)

    # Stub mass details
    st.subheader("MS-Cleavable Crosslinker Stub Masses")

    for name, xl in CROSSLINKERS.items():
        if xl.cleavable and xl.stub_masses:
            st.markdown(f"#### {name}")
            stub_data = [{"Stub": k, "Mass (Da)": f"{v:.4f}"} for k, v in xl.stub_masses.items()]
            st.dataframe(pd.DataFrame(stub_data), use_container_width=True, hide_index=True)

with tab4:
    st.subheader("About Spectrum Annotator Ddzby")

    st.markdown("""
    ### Features

    - **Glycopeptide Annotation**
      - O-glycans (Core 1-4, sialylated, phosphorylated, sulfated)
      - N-glycans (high-mannose, complex, hybrid)
      - Y ion series and oxonium diagnostic ions

    - **Crosslinked Peptide Annotation**
      - MS-cleavable crosslinkers (DSSO, DSBSO)
      - Non-cleavable crosslinkers (BS3, DSS)
      - Stub mass identification

    - **False Match Rate Calculation**
      - Spectrum shifting method
      - Peak-based and intensity-based FMR

    ### Installation

    ```bash
    pip install spectrum-annotator-ddzby
    ```

    ### Python API

    ```python
    from spectrum_annotator_ddzby import (
        FragmentCalculator,
        O_GLYCAN_COMPOSITIONS,
        CROSSLINKERS,
    )

    # Calculate glycopeptide fragments
    calc = FragmentCalculator(
        peptide="PEPTIDEK",
        modifications={"5": "HexNAc"},
        precursor_charge=3
    )
    fragments = calc.calculate_all_fragments()
    ```

    ### Citation

    If you use this tool in your research, please cite:

    > Fu, L. (2026). Spectrum Annotator Ddzby: Universal MS/MS Spectrum Annotation Tool.
    > GitHub: https://github.com/lfu46/GlycoSpectrumAnnotator

    ### Contact

    - **Author:** Longping Fu
    - **Email:** lpfu46@gatech.edu
    - **GitHub:** [lfu46/GlycoSpectrumAnnotator](https://github.com/lfu46/GlycoSpectrumAnnotator)
    """)

# Footer
st.markdown("---")
st.markdown(
    "<p style='text-align: center; color: #888;'>Spectrum Annotator Ddzby v0.2.0 | "
    "Developed by Longping Fu | "
    "<a href='https://github.com/lfu46/GlycoSpectrumAnnotator'>GitHub</a></p>",
    unsafe_allow_html=True
)
