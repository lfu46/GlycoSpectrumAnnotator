#!/usr/bin/env python3
"""
Spectrum Annotator Ddzby - Streamlit Web Application

Interactive web interface for MS/MS spectrum annotation of glycopeptides
and crosslinked peptides.

Inspired by IPSA (Interactive Peptide Spectrum Annotator) from Coon Lab.

Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io

# Try to import PIL for TIFF export
try:
    from PIL import Image
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False

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

# Initialize session state for theme
if 'dark_mode' not in st.session_state:
    st.session_state.dark_mode = False

# IPSA-inspired color scheme
IPSA_COLORS = {
    'header_bg': '#3c4b64',
    'panel_bg': '#3c4b64',
    'b_ion': '#0d75bc',
    'y_ion': '#be202d',
    'c_ion': '#07a14a',
    'z_ion': '#f79420',
    'Y_ion': '#8491B4',
    'oxonium': '#7E6148',
    'precursor': '#666666',
    'unassigned': '#a6a6a6',
}

# Custom CSS based on theme
def get_custom_css(dark_mode):
    # Common styles
    common_css = """
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

        .author-logo {
            text-decoration: none;
        }
        .author-logo:hover {
            opacity: 0.85;
        }
        .author-name {
            font-family: 'Inter', sans-serif;
            font-weight: 600;
            font-size: 1.1rem;
            background: linear-gradient(135deg, #A51C30 0%, #C9102F 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }
        .sub-header {
            font-size: 1.25rem;
            text-align: center;
            margin-bottom: 1.5rem;
        }

        /* Larger tab font */
        .stTabs [data-baseweb="tab"] {
            font-size: 1.05rem !important;
        }

        /* Sidebar styling - compact with larger font */
        div[data-testid="stSidebar"] .stMarkdown p {
            font-size: 1rem !important;
            margin-bottom: 0.3rem !important;
        }
        div[data-testid="stSidebar"] .stMarkdown h1,
        div[data-testid="stSidebar"] .stMarkdown h2,
        div[data-testid="stSidebar"] .stMarkdown h3 {
            font-size: 1.1rem !important;
            margin-top: 0.8rem !important;
            margin-bottom: 0.4rem !important;
        }
        div[data-testid="stSidebar"] label {
            font-size: 0.95rem !important;
        }
        div[data-testid="stSidebar"] .stRadio > div {
            gap: 0.5rem !important;
        }
        div[data-testid="stSidebar"] .stCheckbox {
            padding: 0.1rem 0 !important;
        }
        div[data-testid="stSidebar"] [data-testid="stVerticalBlock"] > div {
            padding-top: 0.2rem !important;
            padding-bottom: 0.2rem !important;
        }
    """

    if dark_mode:
        return f"""
        <style>
            {common_css}
            .header-container {{
                display: flex;
                justify-content: space-between;
                align-items: center;
                background: linear-gradient(90deg, #3c4b64 0%, #2c3b54 100%);
                padding: 0.8rem 1.5rem;
                margin: -1rem -1rem 1rem -1rem;
                border-radius: 0 0 10px 10px;
            }}
            .main-title {{
                font-size: 1.8rem;
                font-weight: bold;
                color: white;
                margin: 0;
            }}
            .author-name {{
                font-family: 'Inter', sans-serif;
                font-weight: 600;
                font-size: 0.95rem;
                color: white;
                line-height: 1.2;
            }}
            .sub-header {{
                color: #aaa;
            }}
            .stTabs [data-baseweb="tab-list"] {{
                gap: 2rem;
                background-color: #2c3b54;
                padding: 0.5rem;
                border-radius: 5px;
            }}
            .stTabs [data-baseweb="tab"] {{
                color: white;
            }}
            div[data-testid="stSidebar"] {{
                background-color: #2c3b54;
            }}
            div[data-testid="stSidebar"] .stMarkdown,
            div[data-testid="stSidebar"] label,
            div[data-testid="stSidebar"] .stSelectbox label,
            div[data-testid="stSidebar"] .stSlider label {{
                color: white !important;
            }}
            .stApp {{
                background-color: #1a1a2e;
            }}
            .stApp > header {{
                background-color: transparent;
            }}
            /* Dark mode text colors */
            .stApp, .stApp p, .stApp span, .stApp label, .stApp div {{
                color: #e0e0e0;
            }}
            .stMarkdown, .stMarkdown p, .stMarkdown h1, .stMarkdown h2, .stMarkdown h3, .stMarkdown h4 {{
                color: #e0e0e0 !important;
            }}
            .stTextInput label, .stNumberInput label, .stSelectbox label, .stTextArea label {{
                color: #e0e0e0 !important;
            }}
            .stMetric label, .stMetric [data-testid="stMetricValue"] {{
                color: #e0e0e0 !important;
            }}
            .stDataFrame {{
                color: #e0e0e0;
            }}
        </style>
        """
    else:
        return f"""
        <style>
            {common_css}
            .header-container {{
                display: flex;
                justify-content: space-between;
                align-items: center;
                background: #f8f9fa;
                padding: 0.8rem 1.5rem;
                margin: -1rem -1rem 1rem -1rem;
                border-bottom: 1px solid #e0e0e0;
            }}
            .main-title {{
                font-size: 1.8rem;
                font-weight: bold;
                color: #2c3e50;
                margin: 0;
            }}
            .sub-header {{
                color: #555;
            }}

            /* Main app background */
            .stApp {{
                background-color: #ffffff;
            }}

            /* Tabs */
            .stTabs [data-baseweb="tab-list"] {{
                gap: 1.5rem;
                background-color: #f0f2f6;
                padding: 0.5rem;
                border-radius: 5px;
            }}
            .stTabs [data-baseweb="tab"] {{
                color: #1a1a1a !important;
            }}

            /* Force all main content text to be dark */
            .main .block-container,
            .main .block-container p,
            .main .block-container span,
            .main .block-container label,
            .main .block-container div,
            .main .stMarkdown,
            .main .stMarkdown p,
            .main .stMarkdown h1,
            .main .stMarkdown h2,
            .main .stMarkdown h3,
            .main .stMarkdown h4,
            .main [data-testid="stMarkdownContainer"],
            .main [data-testid="stMarkdownContainer"] p,
            [data-testid="stAppViewBlockContainer"] p,
            [data-testid="stAppViewBlockContainer"] span,
            [data-testid="stAppViewBlockContainer"] label,
            [data-testid="stVerticalBlock"] p,
            [data-testid="stVerticalBlock"] label,
            [data-testid="stVerticalBlock"] span {{
                color: #1a1a1a !important;
            }}

            /* All labels dark */
            .main label,
            .main .stTextInput label,
            .main .stNumberInput label,
            .main .stSelectbox label,
            .main .stRadio label,
            .main .stCheckbox label,
            .main .stTextArea label,
            .stTextInput label,
            .stNumberInput label,
            .stSelectbox label,
            .stRadio label,
            .stCheckbox label {{
                color: #1a1a1a !important;
            }}

            /* Subheaders */
            .main [data-testid="stSubheader"],
            .main h1, .main h2, .main h3, .main h4 {{
                color: #1a1a1a !important;
            }}

            /* Input fields - light backgrounds with dark text */
            .stTextInput > div > div > input,
            .stNumberInput > div > div > input,
            .stTextArea > div > div > textarea {{
                background-color: #ffffff !important;
                color: #1a1a1a !important;
                border: 1px solid #ced4da !important;
            }}
            .stSelectbox > div > div {{
                background-color: #ffffff !important;
                border: 1px solid #ced4da !important;
            }}
            .stSelectbox > div > div > div,
            .stSelectbox [data-baseweb="select"] span {{
                color: #1a1a1a !important;
            }}

            /* Radio buttons text */
            .stRadio > div > label,
            .stRadio [data-baseweb="radio"] + div {{
                color: #1a1a1a !important;
            }}

            /* Metrics */
            .stMetric label,
            .stMetric [data-testid="stMetricValue"],
            .stMetric [data-testid="stMetricLabel"] {{
                color: #1a1a1a !important;
            }}

            /* Info boxes */
            .stAlert {{
                background-color: #e7f3ff;
                color: #1a1a1a !important;
            }}
            .stAlert p {{
                color: #1a1a1a !important;
            }}

            /* Data editor */
            [data-testid="stDataFrame"] {{
                background-color: #ffffff;
            }}
            [data-testid="stDataFrame"] * {{
                color: #1a1a1a !important;
            }}

            /* Column headers */
            [data-testid="column"] > div > div > div {{
                color: #1a1a1a !important;
            }}
        </style>
        """

st.markdown(get_custom_css(st.session_state.dark_mode), unsafe_allow_html=True)

# Header with logo (IPSA-style layout with author branding)
# Website URL - update when your personal website is live
AUTHOR_WEBSITE = "https://longpingfu.com"  # Placeholder URL

st.markdown(f'''
<div class="header-container">
    <div class="main-title">Spectrum Annotator Ddzby</div>
    <a href="{AUTHOR_WEBSITE}" target="_blank" class="author-logo">
        <span class="author-name">Longping Fu</span>
    </a>
</div>
''', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Universal MS/MS Spectrum Annotation for Glycopeptides and Crosslinked Peptides</p>', unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.header("Settings")

    # Theme toggle (IPSA-inspired)
    st.subheader("Appearance")
    col_theme1, col_theme2 = st.columns(2)
    with col_theme1:
        if st.button("Light", use_container_width=True,
                     type="secondary" if st.session_state.dark_mode else "primary"):
            st.session_state.dark_mode = False
            st.rerun()
    with col_theme2:
        if st.button("Dark", use_container_width=True,
                     type="primary" if st.session_state.dark_mode else "secondary"):
            st.session_state.dark_mode = True
            st.rerun()

    st.divider()

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

    # Ion type toggles (IPSA-style)
    st.subheader("Ion Types")

    st.markdown("**N-terminal**")
    col_n1, col_n2, col_n3 = st.columns(3)
    with col_n1:
        ion_a = st.checkbox("a", value=False)
    with col_n2:
        ion_b = st.checkbox("b", value=True)
    with col_n3:
        ion_c = st.checkbox("c", value=frag_method in ["ETD", "EThcD"])

    st.markdown("**C-terminal**")
    col_c1, col_c2, col_c3 = st.columns(3)
    with col_c1:
        ion_x = st.checkbox("x", value=False)
    with col_c2:
        ion_y = st.checkbox("y", value=True)
    with col_c3:
        ion_z = st.checkbox("z", value=frag_method in ["ETD", "EThcD"])

    st.markdown("**Neutral Losses**")
    col_nl1, col_nl2, col_nl3 = st.columns(3)
    with col_nl1:
        nl_h2o = st.checkbox("-H2O", value=False)
    with col_nl2:
        nl_nh3 = st.checkbox("-NH3", value=False)
    with col_nl3:
        nl_co2 = st.checkbox("-CO2", value=False)

    st.divider()

    # Color customization (IPSA-style)
    st.subheader("Colors")
    col_color1, col_color2 = st.columns(2)
    with col_color1:
        precursor_color = st.color_picker("Precursor", IPSA_COLORS['precursor'])
    with col_color2:
        unassigned_color = st.color_picker("Unassigned", IPSA_COLORS['unassigned'])

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

        # Spectrum input (IPSA-style table)
        st.subheader("Enter Spectral Data Here")

        # Initialize default spectrum data
        if 'spectrum_df' not in st.session_state:
            st.session_state.spectrum_df = pd.DataFrame({
                'Mass To Charge': [204.087, 366.140, 500.25, 650.30, 800.40, None, None, None, None, None],
                'Intensity': [1000.0, 800.0, 500.0, 300.0, 200.0, None, None, None, None, None]
            })

        # Editable data table (IPSA-style)
        edited_df = st.data_editor(
            st.session_state.spectrum_df,
            num_rows="dynamic",
            use_container_width=True,
            height=200,
            column_config={
                "Mass To Charge": st.column_config.NumberColumn(
                    "Mass To Charge",
                    help="m/z value",
                    format="%.4f",
                    min_value=0.0,
                ),
                "Intensity": st.column_config.NumberColumn(
                    "Intensity",
                    help="Peak intensity",
                    format="%.1f",
                    min_value=0.0,
                ),
            },
            hide_index=True,
        )

        # Clear table button
        col_clear, col_space = st.columns([1, 3])
        with col_clear:
            if st.button("Clear table", type="secondary"):
                st.session_state.spectrum_df = pd.DataFrame({
                    'Mass To Charge': [None] * 10,
                    'Intensity': [None] * 10
                })
                st.rerun()

    with col2:
        st.subheader("Results")

        if st.button("Annotate Spectrum", type="primary"):
            # Parse spectrum from data editor
            peaks = []
            for _, row in edited_df.iterrows():
                mz = row['Mass To Charge']
                intensity = row['Intensity']
                if pd.notna(mz) and pd.notna(intensity) and mz > 0 and intensity > 0:
                    peaks.append((float(mz), float(intensity)))

            if peaks:
                mz_values = [p[0] for p in peaks]
                intensities = [p[1] for p in peaks]
                max_intensity = max(intensities) if intensities else 1

                # Normalize to relative abundance (%)
                rel_intensities = [i / max_intensity * 100 for i in intensities]

                # Track matched ions and their errors for error plot
                matched_ions = []
                matched_errors = []
                matched_mz = []

                # Determine plot colors based on theme
                if st.session_state.dark_mode:
                    plot_bg = '#1e1e1e'
                    paper_bg = '#1e1e1e'
                    grid_color = '#444'
                    text_color = 'white'
                    peak_color = unassigned_color
                else:
                    plot_bg = 'white'
                    paper_bg = 'white'
                    grid_color = '#e0e0e0'
                    text_color = 'black'
                    peak_color = unassigned_color

                # Create subplot with main spectrum and error plot (IPSA-style)
                fig = make_subplots(
                    rows=2, cols=1,
                    row_heights=[0.75, 0.25],
                    shared_xaxes=True,
                    vertical_spacing=0.08
                )

                # Add peaks to main spectrum
                for mz, rel_int in zip(mz_values, rel_intensities):
                    fig.add_trace(go.Scatter(
                        x=[mz, mz],
                        y=[0, rel_int],
                        mode='lines',
                        line=dict(color=peak_color, width=1),
                        showlegend=False,
                        hovertemplate=f'm/z: {mz:.4f}<br>Rel. Int: {rel_int:.1f}%<extra></extra>'
                    ), row=1, col=1)

                # Add oxonium ion annotations for glycopeptides
                if analysis_type == "Glycopeptide":
                    for ion_name, ion_mz in OXONIUM_IONS_EXTENDED.items():
                        for i, (mz, rel_int) in enumerate(zip(mz_values, rel_intensities)):
                            if tolerance_unit == "ppm":
                                error_ppm = (mz - ion_mz) / ion_mz * 1e6
                                match = abs(error_ppm) < tolerance
                            else:
                                error_ppm = (mz - ion_mz) / ion_mz * 1e6
                                match = abs(mz - ion_mz) < tolerance

                            if match:
                                # Add colored peak
                                fig.add_trace(go.Scatter(
                                    x=[mz, mz],
                                    y=[0, rel_int],
                                    mode='lines',
                                    line=dict(color=IPSA_COLORS['oxonium'], width=2),
                                    showlegend=False,
                                ), row=1, col=1)
                                # Add annotation
                                fig.add_annotation(
                                    x=mz, y=rel_int,
                                    text=ion_name,
                                    showarrow=True,
                                    arrowhead=2,
                                    arrowsize=1,
                                    arrowwidth=1,
                                    ax=0, ay=-25,
                                    font=dict(size=9, color=IPSA_COLORS['oxonium']),
                                    row=1, col=1
                                )
                                matched_ions.append(ion_name)
                                matched_errors.append(error_ppm)
                                matched_mz.append(mz)
                                break

                # Add error plot markers (IPSA-style)
                if matched_mz:
                    fig.add_trace(go.Scatter(
                        x=matched_mz,
                        y=matched_errors,
                        mode='markers',
                        marker=dict(color=IPSA_COLORS['oxonium'], size=6),
                        showlegend=False,
                        hovertemplate='m/z: %{x:.4f}<br>Error: %{y:.2f} ppm<extra></extra>'
                    ), row=2, col=1)

                # Add zero line to error plot
                fig.add_hline(y=0, line_dash="dash", line_color="gray", row=2, col=1)
                fig.add_hline(y=tolerance, line_dash="dot", line_color="lightgray", row=2, col=1)
                fig.add_hline(y=-tolerance, line_dash="dot", line_color="lightgray", row=2, col=1)

                # Display peptide sequence at top (IPSA-style)
                peptide_display = "   ".join(list(peptide))
                fig.add_annotation(
                    x=0.5, y=1.15,
                    xref="paper", yref="paper",
                    text=f"<b>{peptide_display}</b>",
                    showarrow=False,
                    font=dict(size=16, color=text_color, family="monospace"),
                    xanchor="center"
                )

                # Add spectrum info (IPSA-style)
                info_text = f"Precursor m/z: {precursor_mz if precursor_mz > 0 else 'N/A'}    "
                info_text += f"Charge: +{precursor_charge}    "
                info_text += f"Matched Ions: {len(matched_ions)}/{len(peaks)}"
                fig.add_annotation(
                    x=0.5, y=1.08,
                    xref="paper", yref="paper",
                    text=info_text,
                    showarrow=False,
                    font=dict(size=11, color=text_color),
                    xanchor="center"
                )

                fig.update_layout(
                    plot_bgcolor=plot_bg,
                    paper_bgcolor=paper_bg,
                    showlegend=False,
                    height=600,
                    margin=dict(t=100, b=50),
                    font=dict(color=text_color),
                )

                fig.update_xaxes(
                    title_text="m/z",
                    gridcolor=grid_color,
                    row=2, col=1
                )
                fig.update_yaxes(
                    title_text="Relative Abundance (%)",
                    gridcolor=grid_color,
                    row=1, col=1
                )
                fig.update_yaxes(
                    title_text="Error (ppm)",
                    gridcolor=grid_color,
                    range=[-tolerance * 1.5, tolerance * 1.5],
                    row=2, col=1
                )

                st.plotly_chart(fig, use_container_width=True)

                # Export options
                st.markdown("**Export Options**")
                col_format, col_dpi = st.columns([2, 1])
                with col_format:
                    export_format = st.selectbox(
                        "Format",
                        ["PDF", "TIFF", "EMF", "SVG", "PNG"],
                        index=0,
                        label_visibility="collapsed"
                    )
                with col_dpi:
                    if export_format == "TIFF":
                        tiff_dpi = st.number_input("DPI", min_value=72, max_value=1200, value=600, step=50)
                    else:
                        tiff_dpi = 600  # default

                # Export buttons
                col_btn1, col_btn2, col_btn3 = st.columns(3)
                with col_btn1:
                    try:
                        if export_format == "PDF":
                            pdf_bytes = fig.to_image(format="pdf", width=1200, height=800)
                            st.download_button(
                                "Download PDF",
                                data=pdf_bytes,
                                file_name="spectrum.pdf",
                                mime="application/pdf"
                            )
                        elif export_format == "TIFF":
                            if PIL_AVAILABLE:
                                # Export as PNG first, then convert to TIFF
                                png_bytes = fig.to_image(format="png", width=int(8*tiff_dpi), height=int(6*tiff_dpi), scale=1)
                                img = Image.open(io.BytesIO(png_bytes))
                                tiff_buffer = io.BytesIO()
                                img.save(tiff_buffer, format="TIFF", dpi=(tiff_dpi, tiff_dpi), compression="tiff_lzw")
                                tiff_buffer.seek(0)
                                st.download_button(
                                    f"Download TIFF ({tiff_dpi} DPI)",
                                    data=tiff_buffer.getvalue(),
                                    file_name="spectrum.tiff",
                                    mime="image/tiff"
                                )
                            else:
                                st.warning("Install Pillow for TIFF export: pip install Pillow")
                        elif export_format == "EMF":
                            # EMF export via SVG (users can convert with Inkscape or other tools)
                            svg_bytes = fig.to_image(format="svg", width=1200, height=800)
                            st.download_button(
                                "Download SVG (convert to EMF)",
                                data=svg_bytes,
                                file_name="spectrum.svg",
                                mime="image/svg+xml",
                                help="EMF: Open SVG in Inkscape and save as EMF"
                            )
                        elif export_format == "SVG":
                            svg_bytes = fig.to_image(format="svg", width=1200, height=800)
                            st.download_button(
                                "Download SVG",
                                data=svg_bytes,
                                file_name="spectrum.svg",
                                mime="image/svg+xml"
                            )
                        elif export_format == "PNG":
                            png_bytes = fig.to_image(format="png", width=1200, height=800, scale=2)
                            st.download_button(
                                "Download PNG",
                                data=png_bytes,
                                file_name="spectrum.png",
                                mime="image/png"
                            )
                    except Exception as e:
                        st.error(f"Export failed. Install kaleido: pip install kaleido")

                with col_btn2:
                    # Export matched data
                    export_data = pd.DataFrame({
                        'm/z': matched_mz,
                        'Ion': matched_ions,
                        'Error (ppm)': matched_errors
                    })
                    st.download_button(
                        "Export Annotations",
                        data=export_data.to_csv(index=False),
                        file_name="matched_ions.csv",
                        mime="text/csv"
                    )
                with col_btn3:
                    # Export peak list
                    peak_data = pd.DataFrame({
                        'm/z': mz_values,
                        'Intensity': intensities,
                        'Relative (%)': rel_intensities
                    })
                    st.download_button(
                        "Export Peak List",
                        data=peak_data.to_csv(index=False),
                        file_name="peak_list.csv",
                        mime="text/csv"
                    )

                # Summary statistics
                st.subheader("Annotation Summary")
                coverage = len(matched_ions) / len(peaks) * 100 if peaks else 0
                col_a, col_b, col_c = st.columns(3)
                with col_a:
                    st.metric("Total Peaks", len(peaks))
                with col_b:
                    st.metric("Annotated", len(matched_ions))
                with col_c:
                    st.metric("Coverage", f"{coverage:.1f}%")
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

    - **IPSA-Inspired Interface**
      - Dark/Light theme toggle
      - Ion type selection
      - Error (ppm) subplot
      - Interactive spectrum visualization

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

    ### Citations

    If you use this tool in your research, please cite:

    > Fu, L. (2026). Spectrum Annotator Ddzby: Universal MS/MS Spectrum Annotation Tool.
    > GitHub: https://github.com/lfu46/GlycoSpectrumAnnotator

    ---

    **Acknowledgments**

    The user interface of this tool is inspired by **IPSA (Interactive Peptide Spectrum Annotator)**
    from the Coon Laboratory at the University of Wisconsin-Madison. If you use similar visualization
    approaches, please also consider citing the IPSA paper:

    > Brademan, D.R., Riley, N.M., Kwiecien, N.W., and Coon, J.J. (2019).
    > **Interactive Peptide Spectral Annotator: A Versatile Web-based Tool for Proteomic Applications.**
    > *Molecular & Cellular Proteomics*, 18(8 Suppl 1), S193-S201.
    > DOI: [10.1074/mcp.TIR118.001209](https://doi.org/10.1074/mcp.TIR118.001209)

    IPSA is available at: [https://coonlabs.com/ipsa/](https://coonlabs.com/ipsa/)

    ---

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
