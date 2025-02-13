import streamlit as st
import streamlit.components.v1 as components
from langchain_core.messages import AIMessage, HumanMessage
from graph import invoke_our_graph
from st_callable_util import (
    get_streamlit_cb,
)  # Utility function to get a Streamlit callback handler with context
from dotenv import load_dotenv
import re
from shared_state import (
    substance_names,
    init_conc,
    input_expr,
    output_expr,
    coefficient,
)

load_dotenv()

st.title("Chemical Assistant")


def extract_table_blocks(content):
    """从内容中提取表格代码块和普通文本，并转换为substance_names和init_conc格式"""
    global substance_names, init_conc

    # 使用正则表达式匹配<table>...</table>标签
    pattern = r"<table.*?>(.*?)</table>"
    table_blocks = re.findall(pattern, content, re.DOTALL)

    # 从表格内容中提取物质名称和浓度
    substance_names = []
    init_conc = {}

    for table in table_blocks:
        # 匹配表格行中的物质名称和浓度
        row_pattern = r"<td.*?>(.*?)</td>\s*<td.*?>(.*?)</td>"
        matches = re.findall(row_pattern, table)

        for substance, conc in matches:
            substance = substance.strip()
            if substance and substance not in ["Substance", "物质"]:
                substance_names.append(substance)
                init_conc[substance] = float(conc)

    # 将原始内容中的表格替换为占位符
    text_content = re.sub(pattern, "[TABLE]", content)
    return text_content, table_blocks, substance_names, init_conc


def render_message(content):
    """渲染消息内容，分别处理表格和文本"""
    text_content, table_blocks, _, _ = extract_table_blocks(content)

    # 首先显示文本内容
    if text_content.strip():
        st.markdown(text_content, unsafe_allow_html=True)

    # 然后显示每个表格，添加基本样式
    for table in table_blocks:
        table_html = f"""
        <div style="overflow-x: auto;">
            <table style="width: 100%; border-collapse: collapse;">
                {table}
            </table>
        </div>
        """
        components.html(table_html, height=200, scrolling=True)


def extract_chemical_expression(content):
    """从chemical_balance节点的响应中提取化学表达式，支持多个表达式"""
    global input_expr, output_expr, coefficient

    # 使用正则表达式匹配多个Python字典格式的表达式
    pattern = r"\{([^}]+)\},\s*\{([^}]+)\},\s*([\d\.\*\-/]+)"
    matches = re.finditer(pattern, content)

    expressions = []
    for match in matches:
        try:
            # 提取并评估三个部分
            input_dict = eval("{" + match.group(1) + "}")
            output_dict = eval("{" + match.group(2) + "}")
            coef = eval(match.group(3))
            expressions.append((input_dict, output_dict, coef))
        except Exception as e:
            st.warning(f"解析表达式时出错: {str(e)}")
            continue

    if expressions:
        # 保存最后一个表达式到全局变量（保持向后兼容）
        input_expr, output_expr, coefficient = expressions[-1]
        return expressions
    return None


def render_chemical_expression(expressions):
    """渲染多个化学表达式"""
    if not expressions:
        return

    if isinstance(expressions[0], tuple):
        # 多个表达式的情况
        for i, (input_dict, output_dict, coef) in enumerate(expressions, 1):
            st.code(
                f"""# 化学反应表达式 {i}
input_expr = {input_dict}
output_expr = {output_dict}
coefficient = {coef}""",
                language="python",
            )
    else:
        # 单个表达式的情况（向后兼容）
        input_dict, output_dict, coef = expressions
        st.code(
            f"""# 化学反应表达式
input_expr = {input_dict}
output_expr = {output_dict}
coefficient = {coef}""",
            language="python",
        )


if "messages" not in st.session_state:
    st.session_state["messages"] = [
        AIMessage(
            content="我是帮助你做化学反应的助手，请输入你的化学反应物质和初始化浓度"
        )
    ]

for msg in st.session_state.messages:
    if type(msg) == AIMessage:
        with st.chat_message("assistant"):
            st.write(msg.content)
    if type(msg) == HumanMessage:
        with st.chat_message("user"):
            st.write(msg.content)

if prompt := st.chat_input():
    st.session_state.messages.append(HumanMessage(content=prompt))
    with st.chat_message("user"):
        st.write(prompt)
    with st.chat_message("assistant"):
        msg_placeholder = st.empty()
        st_callback = get_streamlit_cb(msg_placeholder)
        full_response = ""
        response = invoke_our_graph(st.session_state.messages, [st_callback])
        skill = response["next"]
        last_msg = response["messages"][-1].content
        st.session_state.messages.append(AIMessage(content=last_msg))
        if skill == "substance_input":
            render_message(last_msg)
        elif skill == "chemical_balance":
            expressions = extract_chemical_expression(last_msg)
            if expressions:
                render_chemical_expression(expressions)
        elif skill == "calculate":
            st.write(last_msg)
            try:
                st.image("1.png", caption="化学平衡计算结果")
            except Exception as e:
                st.error(f"无法加载图片：{str(e)}")
