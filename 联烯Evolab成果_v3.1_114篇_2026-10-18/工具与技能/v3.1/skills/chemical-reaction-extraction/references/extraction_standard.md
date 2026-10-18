# 用户指定的化学文献反应提取规范

来源：`化学文献反应提取(1).pdf`（3 页，含蓝色批注）及 `反应分类标准-20261016.pdf`（1 页）。源文件哈希见 `standard_sources.json`。以下是用户要求的提取契约，不是论文中的化学事实。

## 1. 分子独立建档

`molecules` 为独立注册表，`reaction.participants` 只通过 `molecule_id` 关联。`molecule_id` 使用不冲突的稳定内部 ID；未核对的分子优先用文献命名空间内的临时 ID。CAS 不是内部主键，也不保证一个物质只有一个 CAS。相同名称、相同编号或不含立体信息的 SMILES 不足以自动合并不同来源的身份。确认同一分子后可以显式追加多篇文献的 `occurrences`，不把文献归属当成身份属性。

| 文档字段组 | v3 字段与要求 |
|---|---|
| 唯一标识 | molecule_id；不能依赖文献中 1a、2a 等重复编号 |
| 名称与标识 | names.original/english/abbreviations/aliases、cas_numbers；不虚构缺失 CAS；英文展开或统一命名注明推导方式 |
| 文献关联 | occurrences.document_id/compound_label/evidence_ids；同一分子可多处引用 |
| 表示方式与结构 | representation_type、smiles 或 generic_template；具体 SMILES 保留所需立体信息 |
| 立体信息 | stereochemistry.atom_centers/double_bonds/axial/mixture/assignment；包括 R/S、E/Z、R_a/S_a、M/P 等原文标记及是否明确；未指定不是没有手性 |
| 通用结构模板 | fixed_structure（确定原子、键、连接及缩写含义）；variables（变量标签、连接位置、原文定义、允许取值、关联/排除约束）；instance_mapping（具体编号、变量赋值、具体 molecule_id）；reason |
| 联烯片段字段 | 按蓝色批注暂不强制手填。allene_features 可选；若后续结构派生，应逐片段保留 C=C=C 原子定位、取代/轴手性和局部环境及依据 |
| 外部性质 | external_properties 只关联身份结构已核对的具体分子；保留库名、库内 ID、匹配方式/状态、原始指标/值/单位/关系符/范围、p 值与计算依据、靶点或被萃取金属。指标不局限 IC50；也包括化学性能。assay_context/statistics 为可选补充内容（schema 中保留空字段以便扩展），不得编造 |

对于尚不能完整用 SMILES 表达的轴手性、混合物或未展开通式，保留模板和明确原文立体描述。模板可以包含未知变量，但不能用 `*` 冒充完整具体分子的 SMILES。允许描述命名混合物的组成未知，不凭常见溶剂名制造唯一组分。

## 2. 反应按具体实验记录

| 文档字段组 | v3 字段与要求 |
|---|---|
| 标识、来源、通式关系 | reaction_id、document_id、record_type、location、general_scheme_ids；具体实验与通式分开 |
| 实验类别 | experiment_categories：optimization、substrate_scope、control、scale_up、derivatization、mechanistic_study、other；真实机理研究实验不等于作者提出的机理 |
| 反应类型 | reaction_name 保留原文命名；transformation_type 及 classification 的自动结果标推导，不冒充原文 |
| 参与组分 | participants.participant_id/molecule_id/role/stage；包含反应物、主副产物、催化剂、配体、试剂、添加剂、溶剂和后处理/纯化组分 |
| 投料和规模 | charges.mass/amount/volume/equivalents/mol_percent；equivalence_reference 明确当量/mol%分母；reaction_scale 保存参照分子及其量 |
| 溶剂 | solvent_conditions.composition/mixture_ratio/total_volume/concentration/concentration_basis；试剂载体与主溶剂分开，投加溶剂量不自动变成混合液实测总体积 |
| 基本条件 | steps 按顺序、阶段保存温度、时间、气氛、压力；不把升温/后处理当成同一次恒温保持 |
| 特殊条件和装置 | special_conditions.light_source/wavelength/power/current/voltage/electrodes/flow_parameters/reactor |
| 实验操作 | procedure.original_text/general_procedure_reference/modifications/addition_order/pretreatment；有 SI 才能读取 SI，不将“通常在 SI 中”解释成已获得 SI |
| 后处理、纯化 | steps.stage=workup/purification；记录淬灭、萃取、洗涤、干燥、柱层析、重结晶及原文参数 |
| 产物和其他结果 | outcomes：target_molecule_ids、product_role、value、yield_type、measurement_basis、scope；kind 区分 yield/product_amount/selectivity/conversion/recovery/no_reaction/trace/decomposition/other。保留多步总收率范围，回收率不是产物收率 |

独立谱图峰值、m/z 及碎片丰度通常放入 `source_context`。有真实制备反应的表征证据仍可关联其反应，不能为表征峰额外创造一条合成。

## 3. 统一字段附属信息

每个关键事实必须包含以下字段，不能只写在 notes 中：

| 要求 | 字段 |
|---|---|
| 原始值 | raw_value（允许原始表格单元格及表头组合、结构图位置等对象；不是伪造的原文 SMILES） |
| 标准化值 | normalized_value；量值保留单位、范围、近似及上下界，不强行将 rt、trace、破折号数值化 |
| 获取方式 | acquisition_method：original_extraction、inherited_general_condition、structure_derivation、external_database、inference；不适用可用 not_applicable |
| 证据位置 | evidence_ids，指向本篇或明确关联 SI 的具体 PDF 页/图表区域；结构推导同样要能回到来源 |
| 数据状态 | data_status：confirmed、pending_review、unreported、unreadable、not_applicable、conflict |
| 核查和推导 | verification：unreviewed/source_checked/inferred/missing；derivation 写明计算、结构映射或推断依据。data_status=confirmed 也不代表实验复现 |
| 冲突 | alternatives 保存双方原始值、标准值和各自证据；conflict 的最终 normalized_value 必须为 null。不同表达不必认定作者数值错误 |

机器导入去除自称 `source_checked`、`source_reviewed`、`structure_checked` 的提升；推导仍保留 `inferred`。人工/助手实际看图的复核产物另存，不能覆盖机器原始记录。JSON/schema 通过不等于化学准确或覆盖完整。

## 4. 形成键分类及组织

精确的 41 类编号与标签见 `bond_taxonomy.json`。分类为多标签，一个实验可同时形成多个键；`classification_index.json` 为类别到 reaction_id 的索引，不能复制成多个实验计数。

记录所有观测转化的原子-原子变化，区分新连接、断裂、升键级、降键级。只给新形成连接/新双三键贴形成标签；原有 C=O、保留的 C=C、由双键降级产生的单键不算新的相应成键类别。C-M 按具体金属（如 C-Pd）标记；未列键按 other-A-B，另保留键级。不遗漏隐式氢的净增减；未确认质子来源或原子映射完整性时不宣称全部原子配对完成。

文档“所有的原子-原子形成的键都列出来”为分类覆盖目标。“三个模型继续做反应分类”为可扩展的分类复核安排：保留 model_votes 槽位，但当前未配置、未执行三模型独立分类，不用同一模型重复回答伪装三模型结果；模型一致不提升为原文核查。

## 验收边界

新 Evolab Session 必须产出 v3 记录并执行 `standard merge`；它自动生成旧兼容视图供既有覆盖/审计工具读取。批处理接收器要求兼容记录 ID 和 v3 reaction_id 一一对应，并验证所有字段、身份引用、证据、分类差分、原文哈希及工具资源未被 Session 改动。旧 v2 历史记录不自动迁移或冒充合规。全量处理仍服从暂停标记，单独样例 Session 不恢复队列。
