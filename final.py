import pandas as pd
import numpy as np
import streamlit as st
import io


def calcular_shannon(datos_filtrados, abundancia):
    shannon = 0
    for conteo in datos_filtrados:
        pi = conteo / abundancia
        shannon += -(pi * np.log(pi))
    return shannon

def calcular_pielou(shannon, riqueza):
    if riqueza > 1:
        return shannon / np.log(riqueza)
    return 0

def calcular_Simpson_D(datos_filtrados, abundancia):
    D = 0
    for conteo in datos_filtrados:
        pi = conteo / abundancia
        D += pi**2
    return D

def calcular_jaccard(columna1, columna2):
    col1_num = pd.to_numeric(columna1, errors='coerce').fillna(0)
    col2_num = pd.to_numeric(columna2, errors='coerce').fillna(0)
    especies_p1 = set(col1_num[col1_num > 0].index)
    especies_p2 = set(col2_num[col2_num > 0].index)
    if len(especies_p1 | especies_p2) == 0:
        return 0.0
    return round(len(especies_p1 & especies_p2) / len(especies_p1 | especies_p2), 4)


st.set_page_config(page_title="BioEco Analizador", layout="wide")

st.title("Analizador Automatizado de Ecologia de Comunidades")
st.markdown("---")

with st.sidebar:
    st.markdown("Sube tu matriz de abundancias en formato Excel o CSV.")
    archivo_subido = st.file_uploader("Selecciona tu archivo", type=["xlsx", "xls", "csv"])
    st.markdown("---")

if archivo_subido is not None:
    resultados_puntos = {}
    lista_sinonimos = ["especie", "especies", "taxon", "taxon", "taxones", "familia", "familias",
                       "bicho", "bichos", "organismo", "organismos",
                       "macroinvertebrado", "macroinvertebrados"]
    try:
        if archivo_subido.name.endswith(".csv"):
            df_bruto = pd.read_csv(archivo_subido, header=None)
        else:
            df_bruto = pd.read_excel(archivo_subido, header=None)

        df_strings = df_bruto.astype(str).apply(lambda x: x.str.lower().str.strip())
        idx_fila, idx_col, palabra_encontrada = None, None, None

        for sinonimo in lista_sinonimos:
            f, c = np.where(df_strings == sinonimo)
            if len(f) > 0:
                idx_fila = int(f[0])
                idx_col  = int(c[0])
                palabra_encontrada = str(df_bruto.iloc[idx_fila, idx_col])
                break

        if palabra_encontrada is None:
            st.error("No se reconocio la columna de identificacion. Renombra esa columna a 'Especie' o 'Taxones'.")
            st.stop()

        cabecera = df_bruto.iloc[idx_fila].tolist()
        df = df_bruto.iloc[idx_fila + 1:].copy()
        df.columns = cabecera

        columnas_validas = [
            col for col in df.columns
            if pd.notna(col) and str(col) != 'nan' and not str(col).startswith('Unnamed')
        ]
        df = df[columnas_validas].copy()

        df.rename(columns={cabecera[idx_col]: palabra_encontrada}, inplace=True)
        df = df[df[palabra_encontrada].notna()]
        df.set_index(palabra_encontrada, inplace=True)

        df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

        st.sidebar.success("Analisis completado con exito")

        for punto in df.columns:
            datos_filtrados = df[punto][df[punto] > 0]
            if len(datos_filtrados) == 0:
                resultados_puntos[punto] = {
                    "Abundancia total (N)": 0, "Riqueza de especies (S)": 0,
                    "Indice de Shannon (H')": 0.0, "Indice de Pielou (J')": 0.0,
                    "Indice de Simpson (D)": 0.0, "Indice de Simpson (1-D)": 0.0,
                    "Indice de Simpson (1/D)": 0.0
                }
                continue
            abundancia = datos_filtrados.sum()
            riqueza    = len(datos_filtrados)
            H = calcular_shannon(datos_filtrados, abundancia)
            J = calcular_pielou(H, riqueza)
            D = calcular_Simpson_D(datos_filtrados, abundancia)
            resultados_puntos[punto] = {
                "Abundancia total (N)":    int(abundancia),
                "Riqueza de especies (S)": int(riqueza),
                "Indice de Shannon (H')":  round(H, 4),
                "Indice de Pielou (J')":   round(J, 4),
                "Indice de Simpson (D)":   round(D, 4),
                "Indice de Simpson (1-D)": round(1 - D, 4),
                "Indice de Simpson (1/D)": round(1 / D, 4) if D != 0 else 0.0
            }

        df_resultados = pd.DataFrame(resultados_puntos)
        df.index.name = None
        df_resultados.index.name = None
        df_final = pd.concat([df, df_resultados], axis=0)
        df_final.index.name = palabra_encontrada

        df_jaccard = pd.DataFrame(index=df.columns, columns=df.columns)
        for p1 in df.columns:
            for p2 in df.columns:
                df_jaccard.loc[p1, p2] = calcular_jaccard(df[p1], df[p2])

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Puntos Analizados", len(df.columns))
        with col2:
            st.metric("Riqueza Total de Taxones", len(df.index))
        with col3:
            st.metric("Taxon Dominante Global", str(df.sum(axis=1).idxmax()))

        st.markdown("---")

        buffer_excel = io.BytesIO()
        with pd.ExcelWriter(buffer_excel, engine='openpyxl') as writer:
            df_final.to_excel(writer, sheet_name="Biodiversidad")
            df_jaccard.to_excel(writer, sheet_name="Similitud_Jaccard")
        buffer_excel.seek(0)

        with st.sidebar:
            st.write("---")
            st.write("### Guardar Resultados")
            st.download_button(
                label="Descargar Excel con Indices",
                data=buffer_excel,
                file_name="resultados_biodiversidad.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

        tab_bio, tab_jac = st.tabs(["Biodiversidad e Indices", "Similitud de Jaccard"])
        with tab_bio:
            st.subheader("Tabla Completa de Datos Calculados")
            st.dataframe(df_final, use_container_width=True)
        with tab_jac:
            st.subheader("Matriz Cuadrada de Jaccard (Presencia/Ausencia)")
            st.dataframe(
                df_jaccard.astype(float).style.background_gradient(cmap="RdYlGn", vmin=0, vmax=1, axis=None),
                use_container_width=True
            )

    except Exception as e:
        st.error(f"Error inesperado: {e}")
        import traceback
        st.code(traceback.format_exc())

else:
    st.info("Bienvenido. Sube tu archivo Excel en el panel izquierdo.")
