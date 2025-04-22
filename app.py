import streamlit as st
from Bio import Entrez
import primer3
import matplotlib.pyplot as plt
import io
import pandas as pd
import plotly.graph_objects as go

# --- 페이지 설정 ---
st.set_page_config(
    page_title="모던 qPCR Primer 디자이너",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- CSS 스타일 적용 ---
st.markdown("""
<style>
    /* 전체 폰트 및 색상 */
    * {font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;}
    
    /* 제목 스타일링 */
    h1 {
        color: #2E86AB;
        font-weight: 700;
        margin-bottom: 1.5rem;
    }
    h2 {
        color: #2E86AB;
        font-weight: 600;
        margin-top: 2rem;
    }
    
    /* 버튼 스타일링 */
    .stButton>button {
        background-color: #2E86AB;
        color: white;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        font-weight: 600;
        border: none;
        width: 100%;
        transition: all 0.3s;
    }
    .stButton>button:hover {
        background-color: #1A6E8E;
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
    }
    
    /* 카드 스타일링 */
    .card {
        background-color: white;
        border-radius: 10px;
        padding: 1.5rem;
        box-shadow: 0 4px 12px rgba(0,0,0,0.05);
        margin-bottom: 1rem;
    }
    
    /* 하이라이트 텍스트 */
    .highlight {
        background-color: #F0F7FA;
        color: #2E86AB;
        padding: 0.2rem 0.5rem;
        border-radius: 4px;
        font-weight: 500;
    }
    
    /* 익스팬더 */
    .streamlit-expanderHeader {
        background-color: #F0F7FA;
        border-radius: 8px;
    }
    
    /* 다운로드 버튼 */
    .download-button {
        background-color: #13C366;
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 8px;
        text-decoration: none;
        font-weight: 600;
        display: inline-block;
        margin-top: 1rem;
    }
    
    /* 분리선 */
    hr {
        margin: 2rem 0;
        border-top: 1px solid #eee;
    }
    
    /* 헤더 영역 */
    .header-container {
        display: flex;
        align-items: center;
        background-color: #F0F7FA;
        padding: 1.5rem;
        border-radius: 12px;
        margin-bottom: 2rem;
    }
    
    /* 뱃지 스타일 */
    .badge {
        background-color: #E3F2FD;
        color: #2E86AB;
        padding: 0.3rem 0.7rem;
        border-radius: 20px;
        font-size: 0.8rem;
        font-weight: 500;
        margin-right: 0.5rem;
    }
    
    /* 결과 테이블 */
    .dataframe {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        font-size: 0.9rem;
    }
    
    /* 로더 스타일 */
    .stSpinner > div {
        border-color: #2E86AB transparent transparent;
    }
    
    /* 사이드바 아이콘 */
    .sidebar-icon {
        font-size: 2.5rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 1rem;
    }
    
    /* 진행 단계 스타일 */
    .step-item {
        display: flex;
        align-items: center;
        margin-bottom: 0.5rem;
    }
    
    .step-number {
        background-color: #2E86AB;
        color: white;
        width: 24px;
        height: 24px;
        border-radius: 50%;
        display: flex;
        align-items: center;
        justify-content: center;
        margin-right: 8px;
        font-size: 0.8rem;
        font-weight: bold;
    }
    
    /* Info 박스 스타일 */
    .info-box {
        background-color: #E3F2FD;
        border-left: 4px solid #2E86AB;
        padding: 1rem;
        border-radius: 4px;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ---------- FUNCTIONS ----------
def fetch_mrna_sequence(gene_name: str, email: str) -> str:
    Entrez.email = email
    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{gene_name}[Gene] AND Homo sapiens[Organism] AND mRNA[Filter]",
        retmode="xml"
    )
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        return None
    seq_id = record["IdList"][0]
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()
    return "".join(fasta_data.split("\n")[1:])  # remove FASTA header

def design_primers(seq_id, sequence, **params):
    return primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': sequence
        },
        {
            'PRIMER_OPT_SIZE': params["opt_size"],
            'PRIMER_MIN_SIZE': params["min_size"],
            'PRIMER_MAX_SIZE': params["max_size"],
            'PRIMER_OPT_TM': params["opt_tm"],
            'PRIMER_MIN_TM': params["min_tm"],
            'PRIMER_MAX_TM': params["max_tm"],
            'PRIMER_MIN_GC': params["min_gc"],
            'PRIMER_MAX_GC': params["max_gc"],
            'PRIMER_MAX_POLY_X': params["max_poly_x"],
            'PRIMER_MAX_END_GC': params["max_gc_3_prime"],
            'PRIMER_PRODUCT_SIZE_RANGE': [[params["product_min"], params["product_max"]]],
            'PRIMER_SALT_MONOVALENT': params["salt_monovalent"],
            'PRIMER_SALT_DIVALENT': params["salt_divalent"],
            'PRIMER_DNTPS': params["dntp_conc"],
            'PRIMER_DNA_CONC': params["primer_conc"]
        }
    )

def plot_primer_positions_plotly(primers, sequence_length, top_n=5):
    fig = go.Figure()
    
    # Add horizontal lines for primer pairs
    for i in range(top_n):
        fig.add_trace(go.Scatter(
            x=[0, sequence_length],
            y=[i, i],
            mode='lines',
            line=dict(color='lightgray', width=4),
            showlegend=False
        ))
        
        try:
            fwd_pos = primers[f"PRIMER_LEFT_{i}"][0]
            fwd_len = primers[f"PRIMER_LEFT_{i}"][1]
            rev_pos = primers[f"PRIMER_RIGHT_{i}"][0]
            rev_len = primers[f"PRIMER_RIGHT_{i}"][1]
            
            # Forward primer
            fig.add_trace(go.Scatter(
                x=[fwd_pos, fwd_pos + fwd_len],
                y=[i, i],
                mode='lines',
                line=dict(color='#2E86AB', width=6),
                name='Forward Primer' if i == 0 else None,
                showlegend=i == 0
            ))
            
            # Reverse primer
            fig.add_trace(go.Scatter(
                x=[rev_pos - rev_len, rev_pos],
                y=[i, i],
                mode='lines',
                line=dict(color='#D31867', width=6),
                name='Reverse Primer' if i == 0 else None,
                showlegend=i == 0
            ))
            
            # Add labels
            fig.add_annotation(
                x=fwd_pos,
                y=i + 0.2,
                text=f"F{i+1}",
                showarrow=False,
                font=dict(color='#2E86AB', size=12)
            )
            
            fig.add_annotation(
                x=rev_pos,
                y=i + 0.2,
                text=f"R{i+1}",
                showarrow=False,
                font=dict(color='#D31867', size=12)
            )
            
        except KeyError:
            break
    
    fig.update_layout(
        title="Primer 위치 시각화",
        xaxis_title="mRNA 위치 (bp)",
        yaxis_title="Primer 쌍",
        template="simple_white",
        height=400,
        margin=dict(l=50, r=50, t=50, b=50),
        yaxis=dict(
            tickvals=list(range(top_n)),
            ticktext=[f"Pair {i+1}" for i in range(top_n)],
            autorange="reversed"
        ),
        xaxis=dict(range=[0, sequence_length]),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    return fig

def get_primer_df(primers, top_n=5):
    rows = []
    for i in range(top_n):
        try:
            rows.append({
                "Pair": f"{i+1}",
                "Forward_Seq": primers[f'PRIMER_LEFT_{i}_SEQUENCE'],
                "Reverse_Seq": primers[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                "Forward_Tm": round(primers[f'PRIMER_LEFT_{i}_TM'], 2),
                "Reverse_Tm": round(primers[f'PRIMER_RIGHT_{i}_TM'], 2),
                "Forward_GC%": round(primers[f'PRIMER_LEFT_{i}_GC_PERCENT'], 2),
                "Reverse_GC%": round(primers[f'PRIMER_RIGHT_{i}_GC_PERCENT'], 2),
                "Amplicon_Size": primers[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
            })
        except KeyError:
            break
    return pd.DataFrame(rows)

# --- 사이드바 ---
with st.sidebar:
    st.markdown("### 🧬 qPCR Primer 설계기")
    st.markdown('<div class="sidebar-icon">🧪</div>', unsafe_allow_html=True)
    
    st.markdown("---")
    st.markdown("### 📋 설계 진행 단계")
    
    # 사이드바 진행 단계를 스타일링된 방식으로 표시
    steps = [
        "유전자 선택",
        "설계 조건 설정",
        "Primer 생성",
        "결과 분석 및 다운로드"
    ]
    
    for i, step in enumerate(steps):
        st.markdown(f"""
        <div class="step-item">
            <div class="step-number">{i+1}</div>
            <div>{step}</div>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("---")
    st.markdown("### 📚 참고 자료")
    st.markdown("[NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)")
    st.markdown("[qPCR 프라이머 설계 가이드](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3525423/)")
    
    st.markdown("---")
    st.caption("© 2025 qPCR Primer Designer")

# ---------- 메인 UI ----------
# 헤더 섹션
st.markdown("""
<div class="header-container">
    <div style="flex: 3">
        <h1>🧬 qPCR Primer 설계기</h1>
        <p>인간 유전자의 qPCR 실험을 위한 최적의 프라이머를 빠르게 설계하세요.</p>
        <div>
            <span class="badge">NCBI</span>
            <span class="badge">Primer3</span>
            <span class="badge">Real-time PCR</span>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

# 검색 섹션
st.markdown('<div class="card">', unsafe_allow_html=True)
col1, col2 = st.columns([3, 1])
with col1:
    gene_name = st.text_input("유전자 이름 또는 심볼 입력", value="TP53", placeholder="예: TP53, BRCA1, EGFR...")
with col2:
    organism = st.selectbox("생물종", ["Homo sapiens", "Mus musculus", "Rattus norvegicus"])
email = "your_email@example.com"  # 사용자 이메일 설정
st.markdown('</div>', unsafe_allow_html=True)

# 탭 생성
tab1, tab2, tab3 = st.tabs(["🔍 프라이머 설계 조건", "⚙️ 고급 설정", "📊 결과 및 시각화"])

# 탭 1: 프라이머 설계 조건
with tab1:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("##### 📏 Primer 길이 설정")
        min_size = st.slider("최소 길이 (bp)", 15, 30, 18, 1)
        opt_size = st.slider("최적 길이 (bp)", 15, 30, 20, 1)
        max_size = st.slider("최대 길이 (bp)", 20, 35, 25, 1)
    
    with col2:
        st.markdown("##### 🌡️ Tm 설정")
        min_tm = st.slider("최소 Tm (°C)", 50.0, 65.0, 58.0, 0.5)
        opt_tm = st.slider("최적 Tm (°C)", 55.0, 70.0, 60.0, 0.5)
        max_tm = st.slider("최대 Tm (°C)", 60.0, 75.0, 62.0, 0.5)
    
    st.markdown("##### 🧪 GC 함량 설정")
    min_gc, max_gc = st.slider("GC 함량 범위 (%)", 30.0, 80.0, (40.0, 60.0), 1.0)
    
    st.markdown("##### 🧬 Amplicon 설정")
    product_min, product_max = st.slider("증폭 산물 길이 범위 (bp)", 50, 500, (80, 150), 5)
    
    st.markdown("</div>", unsafe_allow_html=True)

# 탭 2: 고급 설정
with tab2:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("##### 🔬 특이 조건 설정")
        max_poly_x = st.slider("최대 연속 동일 염기 (PolyX)", 2, 10, 4, 1)
        max_gc_3_prime = st.slider("3' 말단 G/C 최대 개수 (5bp 내)", 0, 5, 3, 1)
    
    with col2:
        st.markdown("##### 🧪 반응 조건")
        salt_monovalent = st.slider("Monovalent Salt (mM)", 0.0, 100.0, 50.0, 1.0)
        salt_divalent = st.slider("Divalent Salt (mM)", 0.0, 10.0, 1.5, 0.1)
        dntp_conc = st.slider("dNTP 농도 (mM)", 0.0, 1.0, 0.2, 0.05)
        primer_conc = st.slider("Primer DNA 농도 (nM)", 10.0, 1000.0, 500.0, 10.0)
    
    st.markdown("##### 🔢 출력 설정")
    top_n = st.slider("출력할 Primer 쌍 개수", 1, 10, 5, 1)
    
    st.markdown("</div>", unsafe_allow_html=True)

# 버튼 섹션
col1, col2, col3 = st.columns([1, 2, 1])
with col2:
    start_button = st.button("🎯 Primer 설계 시작", use_container_width=True)

# 탭 3: 결과 및 시각화 (초기에는 비어있음)
with tab3:
    if not start_button:
        st.markdown("""
        <div class="info-box">
            <h4>👇 Primer 설계를 시작하려면 버튼을 클릭하세요</h4>
            <p>설계 조건을 모두 설정한 후, 아래 '🎯 Primer 설계 시작' 버튼을 클릭하면 이곳에 결과가 표시됩니다.</p>
        </div>
        """, unsafe_allow_html=True)

# 결과 처리
if start_button:
    # 포커스 자동으로 결과 탭으로 이동
    
    
    with st.spinner("🔍 NCBI에서 mRNA 서열을 검색하는 중..."):
        # 로티 애니메이션 대신 진행 상황 표시
        progress_bar = st.progress(0)
        for percent_complete in range(100):
            progress_bar.progress(percent_complete + 1)
            if percent_complete == 50:  # 50% 진행 시점에서 서열 가져오기
                sequence = fetch_mrna_sequence(gene_name, email)
    
    if not sequence:
        st.error("❌ 해당 유전자의 mRNA 서열을 찾을 수 없습니다. 유전자 이름을 확인해주세요.")
    else:
        with tab3:
            # 성공 메시지
            st.success(f"✅ mRNA 서열을 성공적으로 가져왔습니다! ({len(sequence)}bp)")
            
            # Primer 설계
            with st.spinner("🧬 Primer3로 최적의 프라이머 설계 중..."):
                progress_bar = st.progress(0)
                for percent_complete in range(100):
                    progress_bar.progress(percent_complete + 1)
                    if percent_complete == 70:  # 70% 진행 시점에서 프라이머 설계
                        primers = design_primers(gene_name, sequence[:2000],
                            min_size=min_size, opt_size=opt_size, max_size=max_size,
                            min_tm=min_tm, opt_tm=opt_tm, max_tm=max_tm,
                            min_gc=min_gc, max_gc=max_gc,
                            max_poly_x=max_poly_x, max_gc_3_prime=max_gc_3_prime,
                            product_min=product_min, product_max=product_max,
                            salt_monovalent=salt_monovalent, salt_divalent=salt_divalent,
                            dntp_conc=dntp_conc, primer_conc=primer_conc
                        )
            
            # 결과 출력
            st.markdown("## 🎯 분석 결과")
            
            # 프라이머 정보 테이블
            st.markdown("### 📋 프라이머 목록")
            df = get_primer_df(primers, top_n=top_n)
            
            if len(df) == 0:
                st.warning("❗ 설정된 조건에 맞는 프라이머를 찾을 수 없습니다. 조건을 완화해보세요.")
            else:
                st.dataframe(df, use_container_width=True)
                
                # 카드 형태로 추천 프라이머 표시
                st.markdown("### 💎 추천 프라이머 상세 정보")
                
                for i in range(min(len(df), top_n)):
                    with st.expander(f"🔹 Primer Pair {i+1} - Amplicon Size: {df.iloc[i]['Amplicon_Size']}bp"):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown("**Forward Primer:**")
                            st.code(df.iloc[i]['Forward_Seq'])
                            st.markdown(f"""
                            - 길이: {len(df.iloc[i]['Forward_Seq'])}bp
                            - Tm: {df.iloc[i]['Forward_Tm']}°C
                            - GC 함량: {df.iloc[i]['Forward_GC%']}%
                            """)
                        
                        with col2:
                            st.markdown("**Reverse Primer:**")
                            st.code(df.iloc[i]['Reverse_Seq'])
                            st.markdown(f"""
                            - 길이: {len(df.iloc[i]['Reverse_Seq'])}bp
                            - Tm: {df.iloc[i]['Reverse_Tm']}°C
                            - GC 함량: {df.iloc[i]['Reverse_GC%']}%
                            """)
                
                # 플롯티로 시각화
                st.markdown("### 📊 프라이머 위치 시각화")
                fig = plot_primer_positions_plotly(primers, len(sequence[:2000]), top_n=top_n)
                st.plotly_chart(fig, use_container_width=True)
                
                # 서열 정보
                with st.expander("🧬 mRNA 서열 정보"):
                    st.text_area("FASTA 서열", sequence[:100] + "..." + sequence[-100:], height=150)
                    st.markdown(f"전체 길이: **{len(sequence)}bp**")
                
                # 다운로드 섹션
                st.markdown("### 📥 결과 다운로드")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    csv = df.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        label="📋 프라이머 결과 CSV 다운로드",
                        data=csv,
                        file_name=f"{gene_name}_primers.csv",
                        mime="text/csv"
                    )
                
                with col2:
                    fasta_content = f">{gene_name} mRNA\n{sequence}"
                    st.download_button(
                        label="🧬 mRNA 서열 FASTA 다운로드",
                        data=fasta_content,
                        file_name=f"{gene_name}_mRNA.fasta",
                        mime="text/plain"
                    )
                
                # 추가 링크
                st.markdown("### 🔗 추가 분석 도구")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown("""
                    <a href="https://www.ncbi.nlm.nih.gov/tools/primer-blast/" target="_blank" class="download-button">
                        NCBI Primer-BLAST
                    </a>
                    """, unsafe_allow_html=True)
               
                
# 푸터 섹션
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666;">
    <p>© 2025 qPCR Primer Designer | 생화학 연구 및 분자 생물학을 위한 도구</p>
</div>
""", unsafe_allow_html=True)