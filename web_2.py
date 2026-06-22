import numpy as np
import pandas as pd
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
# CARGA EL EXCEL O CSV, BUSCA LA COLUMNA ADECUADA CON TAXONES, LIMPIA Y DEVUELVE UN ARCHIVO LIMPIO
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

    # Busca la palabra clave para indexar. Pasa el texto a minúsculas y elimina espacios para evitar errores de coincidencia
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
# CALCULA TODOS LOS INDICES ECOLÓGICOS
# Ahora acepta parámetros booleanos para calcular solo los índices seleccionados.
# Shannon se calcula siempre que Pielou esté activo, ya que Pielou depende de él.
# ==========================================
def calcular_matriz_completa(df_base, palabra_encontrada, shannon=True, pielou=True, simpson=True, jaccard=True):
    resultados_puntos = {}
    
    # Si Pielou está activo, Shannon debe calcularse aunque el usuario no lo haya marcado explícitamente
    shannon_necesario = shannon or pielou

    for punto in df_base.columns:
        datos_filtrados = df_base[punto][df_base[punto] > 0]
        
        # Siempre se incluyen abundancia y riqueza como base informativa
        resultado = {
            "Abundancia total (N)": 0,
            "Riqueza de especies (S)": 0,
        }

        if len(datos_filtrados) == 0:
            resultados_puntos[punto] = resultado
            continue

        abundancia = datos_filtrados.sum()
        riqueza = len(datos_filtrados)

        resultado["Abundancia total (N)"] = int(abundancia)
        resultado["Riqueza de especies (S)"] = int(riqueza)

        # Calcula Shannon si está marcado o si lo necesita Pielou
        H = None
        if shannon_necesario:
            H = calcular_shannon(datos_filtrados, abundancia)
            if shannon:
                resultado["Indice de Shannon (H')"] = round(H, 4)

        # Pielou requiere Shannon, por eso H nunca será None aquí si pielou es True
        if pielou and H is not None:
            J = calcular_pielou(H, riqueza)
            resultado["Indice de Pielou (J')"] = round(J, 4)

        if simpson:
            D = calcular_Simpson_D(datos_filtrados, abundancia)
            resultado["Indice de Simpson (D)"]   = round(D, 4)
            resultado["Indice de Simpson (1-D)"] = round(1 - D, 4)
            resultado["Indice de Simpson (1/D)"] = round(1 / D, 4) if D != 0 else np.nan

        resultados_puntos[punto] = resultado

    df_indices = pd.DataFrame(resultados_puntos)

    # MATRIZ JACCARD — solo se calcula si el usuario la ha seleccionado
    df_jaccard = None
    if jaccard:
        df_jaccard = pd.DataFrame(index=df_base.columns, columns=df_base.columns)
        for p1 in df_base.columns:
            for p2 in df_base.columns:
                df_jaccard.loc[p1, p2] = calcular_jaccard(df_base[p1], df_base[p2])

    return df_base, df_indices, df_jaccard


# ==========================================
# INTERFAZ DE LA WEB CON STREAMLIT
# ==========================================

st.set_page_config(page_title="BioEco Analizador", layout="wide")

# ES LA PÁGINA PRINCIPAL
st.title("Analizador Automatizado de Ecología de Comunidades")
st.markdown("---")

# ESTO ES LA BARRA LATERAL
with st.sidebar:
    st.markdown("Sube tu matriz de abundancias en formato Excel o CSV.")
    archivo_subido = st.file_uploader("Selecciona tu archivo", type=["xlsx", "xls", "csv"])
    st.markdown("<hr style='margin: 8px 0'>", unsafe_allow_html=True)

    # ------------------------------------------
    # SELECTOR DE ÍNDICES
    # Solo se muestra si hay un archivo cargado.
    # El checkbox "Seleccionar todos" actúa como maestro: si está marcado,
    # fuerza a True el valor inicial de todos los demás checkboxes.
    # Pielou fuerza Shannon a True y lo deshabilita para evitar inconsistencias.
    # ------------------------------------------
    
    st.markdown("### Índices a calcular")

    todos = st.checkbox("Seleccionar todos", value=True)

    calc_pielou  = st.checkbox("Índice de Pielou (J')",            value=todos, help="Al activarse activaautomáticamente el índice de Shannon.")

    # Shannon se deshabilita y fuerza a True si Pielou está marcado,
    # porque Pielou depende matemáticamente de Shannon
    calc_shannon = st.checkbox(
        "Índice de Shannon (H')",
        value=todos or calc_pielou,
        disabled=calc_pielou,
        help="Se calcula usando el logaritmo en base de e."
    )

    calc_simpson = st.checkbox("Índices de Simpson (D, 1-D, 1/D)", value=todos)
    calc_jaccard = st.checkbox("Similitud de Jaccard",              value=todos)

    st.markdown("<hr style='margin: 8px 0'>", unsafe_allow_html=True)
        
    if archivo_subido is not None:


        # ------------------------------------------
        # BOTÓN DE CALCULAR
        # El cálculo ya no ocurre automáticamente al subir el archivo.
        # El usuario elige sus índices y luego pulsa este botón.
        # ------------------------------------------
        calcular = st.button("Calcular", use_container_width=True)

    st.markdown("Desarrollado por Ander MG · 2026")


# ESTO VUELVE A SER PARTE DE LO PRINCIPAL

if archivo_subido is not None:
    try:
        df_base, palabra_encontrada = cargar_y_limpiar_matriz(archivo_subido)

        if palabra_encontrada is None:
            st.error("No se reconoció la columna de identificación. Renombra esa columna a 'Especie' o 'Taxones'.")
            st.stop()

        # ------------------------------------------
        # El bloque de resultados solo se ejecuta cuando el usuario pulsa "Calcular".
        # Antes de pulsarlo, la web queda en espera mostrando un mensaje informativo.
        # ------------------------------------------
        if not calcular:
            st.info("Archivo cargado correctamente. Selecciona los índices en el panel izquierdo y pulsa **Calcular**.")
            st.stop()

        # Se pasan al cálculo solo los índices que el usuario ha marcado
        df_abundancia, df_indices, df_jaccard = calcular_matriz_completa(
            df_base, palabra_encontrada,
            shannon=calc_shannon,
            pielou=calc_pielou,
            simpson=calc_simpson,
            jaccard=calc_jaccard
        )

        st.success("Análisis completado con éxito")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Puntos Analizados", len(df_base.columns))
        with col2:
            st.metric("Riqueza Total de Taxones", len(df_base.index))
        with col3:
            st.metric("Taxon Dominante Global", str(df_base.sum(axis=1).idxmax()))

        st.markdown("---")

        # ------------------------------------------
        # GENERACIÓN DEL EXCEL DESCARGABLE
        # Si Jaccard no fue calculado, su hoja no se incluye en el Excel.
        # ------------------------------------------
        buffer_excel = io.BytesIO()
        with pd.ExcelWriter(buffer_excel, engine='openpyxl') as writer:
            df_abundancia.to_excel(writer, sheet_name="Abundancia")
            df_indices.to_excel(writer, sheet_name="Indices_ecologicos")
            if df_jaccard is not None:
                df_jaccard.to_excel(writer, sheet_name="Similitud_Jaccard")
        buffer_excel.seek(0)

        # ------------------------------------------
        # PESTAÑAS DE RESULTADOS
        # La pestaña de Jaccard solo aparece si el usuario la seleccionó.
        # ------------------------------------------
        tabs_disponibles = ["Biodiversidad e Índices"]
        if df_jaccard is not None:
            tabs_disponibles.append("Similitud de Jaccard")

        tabs = st.tabs(tabs_disponibles)

        with tabs[0]:
            st.subheader("Tabla Completa de Datos Calculados")
            st.download_button(
                label="📥 Descargar Excel con Índices",
                data=buffer_excel,
                file_name="resultados_biodiversidad.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            st.write("")
            st.dataframe(df_indices, use_container_width=True)

        # Solo se renderiza la pestaña de Jaccard si existe
        if df_jaccard is not None:
            with tabs[1]:
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
#  python -m streamlit run web_2.py
#