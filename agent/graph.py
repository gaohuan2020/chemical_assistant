from gc import callbacks
from typing import Annotated
from sympy import im
from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END, MessagesState
from langgraph.graph.message import add_messages
from prompt import substance_input_prompt, chemical_balance_prompt
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from typing import Literal
from langgraph.types import Command
import shared_state


class State(TypedDict):
    # Messages have the type "list". The `add_messages` function
    # in the annotation defines how this state key should be updated
    # (in this case, it appends messages to the list, rather than overwriting them)
    messages: Annotated[list, add_messages]


graph_builder = StateGraph(State)
from langchain_openai import ChatOpenAI

llm = ChatOpenAI(
    openai_api_base="https://api.siliconflow.cn/v1/",
    model="deepseek-ai/DeepSeek-V3",
    openai_api_key="sk-grtxiqhiihsrfiftgjjgrxroyqpwkjnmwukjoickhvdarvig",
    max_tokens=4000,
    temperature=0,
    streaming=True,
)
members = ["substance_input", "chemical_balance", "calculate"]
options = members + ["FINISH"]
system_prompt = (
    "You are a supervisor tasked with managing a conversation between the"
    f" following workers: {members}. Given the following user request,"
    " respond with the worker to act next. If user request is not clear, respond with chitchat. "
    " Each worker will perform a task and respond with their results and status. When finished,"
    " respond with FINISH. The output must be JSON format."
    "EXAMPLE Input: "
    " 定义参与反应的物质和初始化浓度"
    "EXAMPLE JSON OUTPUT: "
    "{"
    '    next: "substance_input" '
    "}"
    "EXAMPLE Input: "
    "将浓度改为30"
    "EXAMPLE JSON OUTPUT: "
    "{"
    '    next: "substance_input" '
    "}"
    "EXAMPLE Input: "
    "定义化学平衡水的电离平衡和铵离子的解离平衡"
    "EXAMPLE JSON OUTPUT: "
    "{"
    '    next: "chemical_balance" '
    "}"
    "EXAMPLE Input: "
    "开始计算"
    "EXAMPLE JSON OUTPUT: "
    "{"
    '    next: "calculate" '
    "}"
)


class Router(TypedDict):
    """Worker to route to next. If no workers needed, route to FINISH."""

    next: Literal[*options]


class State(MessagesState):
    next: str


def calculate(substance_names, init_conc):
    import numpy as np
    import matplotlib.pyplot as plt
    from library.chemistry import Species, Equilibrium
    from library.equilibria import EqSystem
    from matplotlib.font_manager import FontManager
    import matplotlib.font_manager as fm

    # 获取系统中文字体
    fm = FontManager()
    font_names = sorted([f.name for f in fm.ttflist])
    # 设置中文字体
    plt.rcParams["font.sans-serif"] = ["SimHei"]  # 用来正常显示中文标签
    plt.rcParams["axes.unicode_minus"] = False  # 用来正常显示负号
    substance_names = ["H+", "OH-", "NH4+", "NH3", "H2O"]
    # 1. 定义参与反应的物质
    subst = {n: Species.from_formula(n) for n in substance_names}
    init_conc = {
        "H+": 1e-7,  # 0.0000001 M
        "OH-": 1e-7,  # 0.0000001 M
        "NH4+": 1e-7,  # 0.0000001 M
        "NH3": 1.0,  # 1 M
        "H2O": 55.5,  # 55.5 M (纯水的浓度)
    }

    # 3. 定义化学平衡
    # 水的电离平衡：H2O ⇌ H+ + OH-
    w_autop = Equilibrium(
        {"H2O": 1}, {"H+": 1, "OH-": 1}, 10**-14 / 55.5  # 反应物  # 生成物  # 平衡常数
    )

    # 铵离子的解离平衡：NH4+ ⇌ H+ + NH3
    NH4p_pr = Equilibrium(
        {"NH4+": 1}, {"H+": 1, "NH3": 1}, 10**-9.26  # 反应物  # 生成物  # 平衡常数
    )

    equilibria = (w_autop, NH4p_pr)

    # 4. 创建平衡系统
    eqsys = EqSystem(equilibria, subst)

    # 5. 计算不同H+浓度下的平衡
    nc = 60  # 取60个点
    Hp_0 = np.logspace(-3, 0, nc)  # H+浓度从0.001到1

    # 6. 绘图
    plt.figure(figsize=(10, 6))
    plt.grid(True)

    # 计算并绘制结果
    results = eqsys.roots(init_conc, Hp_0, "H+")
    x, info, sane = results

    # 设置最小浓度阈值，避免对数刻度的问题
    min_conc = 1e-20  # 设置一个极小的浓度值作为下限
    x = np.maximum(x, min_conc)  # 将所有小于min_conc的值替换为min_conc

    # 在绘图之前添加调试信息
    print("数据范围：")
    for i, name in enumerate(substance_names):
        print(f"{name}: {x[:, i].min():.2e} to {x[:, i].max():.2e}")

    # 在添加pH刻度之前
    print("pH范围：", min(-np.log10(Hp_0)), "到", max(-np.log10(Hp_0)))

    # 为每个物质画一条线
    for i, name in enumerate(substance_names):
        plt.plot(Hp_0, x[:, i], label=name, linewidth=2)

    # 设置图表格式
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("H+ 浓度 (M)", fontsize=12)
    plt.ylabel("物质浓度 (M)", fontsize=12)
    plt.title("氨水体系中各物质浓度随pH变化的关系", fontsize=14)
    plt.legend(fontsize=10, loc="best")
    plt.grid(True, which="both", ls="-", alpha=0.2)

    # 添加pH刻度
    ax2 = plt.gca().twiny()
    pH_values = [-np.log10(h) for h in Hp_0]
    ax2.set_xlim(plt.gca().get_xlim())

    plt.tight_layout()
    plt.show()
    plt.savefig("1.png")


def supervisor_node(state: State) -> Command[Literal[*members, "__end__"]]:
    print("====================supervisor_node=====================")
    messages = [
        {"role": "system", "content": system_prompt},
    ] + state["messages"]
    response = llm.with_structured_output(Router).invoke(messages)
    goto = response["next"]
    if goto == "FINISH":
        goto = END

    return Command(goto=goto, update={"next": goto})


def substance_input_bot(state: State):
    print("====================substance_node=====================")
    result = llm.invoke(
        [SystemMessage(content=substance_input_prompt)] + state["messages"]
    )
    print("====================Finish substance_node====================")
    return Command(
        update={
            "messages": [HumanMessage(content=result.content, name="substance_input")]
        },
        goto=END,
    )


def chemical_balance_bot(state: State):
    print("====================chemical_balance_node=====================")
    result = llm.invoke(
        [SystemMessage(content=chemical_balance_prompt)] + state["messages"]
    )
    print("====================Finish chemical_balance_node====================")
    return Command(
        update={
            "messages": [HumanMessage(content=result.content, name="chemical_balance")]
        },
        goto=END,
    )


def calculate_bot(state: State):
    print("====================calculate_bot=====================")
    calculate(shared_state.substance_names, shared_state.init_conc)
    print("====================Finish calculate_bot====================")
    return Command(
        update={
            "messages": [HumanMessage(content="计算完成", name="chemical_balance")]
        },
        goto=END,
    )


# The first argument is the unique node name
# The second argument is the function or object that will be called whenever
# the node is used.
builder = StateGraph(State)
builder.add_edge(START, "supervisor")
builder.add_node("supervisor", supervisor_node)
builder.add_node("chemical_balance", chemical_balance_bot)
builder.add_node("substance_input", substance_input_bot)
builder.add_node("calculate", calculate_bot)
builder.add_edge("substance_input", END)
builder.add_edge("supervisor", END)
graph = builder.compile()


def invoke_our_graph(st_messages, callables):
    if not isinstance(callables, list):
        raise TypeError("callables must be a list")

    # Convert all history messages to tuples of (role, content)
    history_messages = []
    for msg in st_messages:
        if isinstance(msg, HumanMessage):
            history_messages.append(("user", msg.content))
        elif isinstance(msg, AIMessage):
            history_messages.append(("assistant", msg.content))

    response = graph.invoke(
        {"messages": history_messages}, subgraphs=True, config={"callbacks": callables}
    )
    last_response = response[-1]
    return last_response


# for response in graph.stream(
#     {
#         "messages": [
#             (
#                 "user",
#                 '参与反应的物质"H+", "OH-", "NH4+", "NH3", "H2O"，对应初始浓度为"H+": 1e-7, # 0.0000001 M "OH-": 1e-7, # 0.0000001 M "NH4+": 1e-7, # 0.0000001 M "NH3": 1.0, # 1 M "H2O": 55.5, # 55.5 M (纯水的浓度)',
#             )
#         ]
#     },
#     subgraphs=True,
# ):
#     last_response = response[-1]
# print(last_response)
# 定义水的电离平衡和铵离子的解离平衡