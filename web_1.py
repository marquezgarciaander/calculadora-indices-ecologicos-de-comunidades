import pandas as pd
import numpy as np
import streamlit as st
import io


# ==========================================
# DEFINIR LOS CÁLCULOS
# ==========================================

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


# ==========================================
# CARGA EL EXCEL O SCV, BUSCA LA COLUMNA ADECUADA CON TAXONES, LIMPIA Y DEVUELVE UN ARCHIVO LIMPIO
# ==========================================
def cargar_y_limpiar_matriz(archivo_subido):
    lista_sinonimos = [
        "especie", "especies", "taxon", "taxones", "familia", "familias",
        "bicho", "bichos", "organismo", "organismos",
        "macroinvertebrado", "macroinvertebrados"
    ]
    
    # Lee el archivo según su extensión
    if archivo_subido.name.endswith(".csv"):
        df_bruto = pd.read_csv(archivo_subido, header=None)
    else:
        df_bruto = pd.read_excel(archivo_subido, header=None)

    # Busca la palabra clave para indexar. Pasa el texto a minusculas y elimina espacios para evitar errores de coincidencia
    df_strings = df_bruto.astype(str).apply(lambda x: x.str.lower().str.strip())
    idx_fila, idx_col, palabra_encontrada = None, None, None

    for sinonimo in lista_sinonimos:
        f, c = np.where(df_strings == sinonimo)
        if len(f) > 0:
            idx_fila = int(f[0])
            idx_col  = int(c[0])
            palabra_encontrada = str(df_bruto.iloc[idx_fila, idx_col])
            break

    # Si no encuentra nada, devolvemos None para manejar el error fuera
    if palabra_encontrada is None:
        return None, None

    # Reestructurar las cabeceras y limpiar datos vacíos
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

    # Convertir datos a enteros numéricos
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    
    return df, palabra_encontrada


# ==========================================
# TCALCULA TODOS LOS INDICES ECOLÓGICOS
# ==========================================
def calcular_matriz_completa(df_base, palabra_encontrada):
    resultados_puntos = {}
    
    for punto in df_base.columns:
        datos_filtrados = df_base[punto][df_base[punto] > 0]
        
        if len(datos_filtrados) == 0:
            resultados_puntos[punto] = {
                "Abundancia total (N)": 0,
                "Riqueza de especies (S)": 0,
                "Indice de Shannon (H')": 0.0,
                "Indice de Pielou (J')": 0.0,
                "Indice de Simpson (D)": 0.0,
                "Indice de Simpson (1-D)": 0.0,
                "Indice de Simpson (1/D)": 0.0
            }
            continue

        abundancia = datos_filtrados.sum()
        riqueza = len(datos_filtrados)

        H = calcular_shannon(datos_filtrados, abundancia)
        J = calcular_pielou(H, riqueza)
        D = calcular_Simpson_D(datos_filtrados, abundancia)

        resultados_puntos[punto] = {
            "Abundancia total (N)": int(abundancia),
            "Riqueza de especies (S)": int(riqueza),
            "Indice de Shannon (H')": round(H, 4),
            "Indice de Pielou (J')": round(J, 4),
            "Indice de Simpson (D)": round(D, 4),
            "Indice de Simpson (1-D)": round(1 - D, 4),
            "Indice de Simpson (1/D)": round(1 / D, 4) if D != 0 else np.nan
        }

    df_indices = pd.DataFrame(resultados_puntos)

    # MATRIZ JACCARD
    df_jaccard = pd.DataFrame(index=df_base.columns, columns=df_base.columns)

    for p1 in df_base.columns:
        for p2 in df_base.columns:
            df_jaccard.loc[p1, p2] = calcular_jaccard(df_base[p1], df_base[p2])

    return df_base, df_indices, df_jaccard


# ==========================================
# INTERFAZ DE LA WEB CON STREAMLIT
# ==========================================

st.set_page_config(page_title="BioEco Analizador", layout="wide")

#ES LA PAGINA PRINCIPAL
st.title("Analizador Automatizado de Ecologia de Comunidades")
st.markdown("---")

#ESTO ES LA BARRA LATERAL
with st.sidebar:
    st.markdown("Sube tu matriz de abundancias en formato Excel o CSV.")
    archivo_subido = st.file_uploader("Selecciona tu archivo", type=["xlsx", "xls", "csv"])
    st.markdown("---")
    st.markdown("Desarrollado por Ander MG · 2026")

#ESTO VUELVE A SER PARTE DE LO PIRNCIPAL
if archivo_subido is not None:
    try:
        # se obtienen los resultados
        df_base, palabra_encontrada = cargar_y_limpiar_matriz(archivo_subido)

        if palabra_encontrada is None:
            st.error("No se reconocio la columna de identificacion. Renombra esa columna a 'Especie' o 'Taxones'.")
            st.stop()
    
        # se obtienen los resultados
        df_abundancia, df_indices, df_jaccard = calcular_matriz_completa(df_base, palabra_encontrada)

        # Muestra el aviso de éxito en el área principal
        st.success("Analisis completado con exito")
   
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Puntos Analizados", len(df_base.columns))
        with col2:
            st.metric("Riqueza Total de Taxones", len(df_base.index))
        with col3:
            st.metric("Taxon Dominante Global", str(df_base.sum(axis=1).idxmax()))

        st.markdown("---")

        buffer_excel = io.BytesIO()
        with pd.ExcelWriter(buffer_excel, engine='openpyxl') as writer:
            df_abundancia.to_excel(writer, sheet_name="Abundancia")
            df_indices.to_excel(writer, sheet_name="Indices_ecologicos")
            df_jaccard.to_excel(writer, sheet_name="Similitud_Jaccard")
        buffer_excel.seek(0)

        tab_bio, tab_jac = st.tabs(["Biodiversidad e Indices", "Similitud de Jaccard"])
        with tab_bio:
            st.subheader("Tabla Completa de Datos Calculados")
            
            # El botón de descarga ahora está integrado en la sección principal de datos
            st.download_button(
                label="📥 Descargar Excel con Indices",
                data=buffer_excel,
                file_name="resultados_biodiversidad.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            st.write("") 
            st.dataframe(df_abundancia, use_container_width=True)
            
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



#
#  python -m streamlit run web_1.py
#