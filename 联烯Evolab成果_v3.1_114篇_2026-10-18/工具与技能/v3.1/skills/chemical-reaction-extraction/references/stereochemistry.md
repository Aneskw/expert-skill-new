# 立体化学与异构体（v3.1）

先查看结构图、楔线/虚线、波浪键、命名、表格脚注和实验/SI的相应出处。原文没指定构型不能解释成无手性、外消旋或单一异构体。原图信息不足时保留候选和具体未解决项。

## 分子身份

`stereochemistry` 的每项仍是带原始值、标准值、证据与状态的 fact：

- `atom_centers`：逐中心记录原文 R/S、相对构型或未知；用 atom_map（若有）或 SMILES 零基 atom_index 定位，注明索引对应的 SMILES。@/@@ 是局部邻居顺序标记，不能直接等同 R/S。
- `double_bonds`：逐双键的 E/Z、cis/trans 和未指定状态，定位两端原子。环的顺反及相对构型存 `assignment`，不强转成双键 E/Z；cis/trans 不无条件等于 Z/E。
- `axial`：联烯、联芳基等轴及原文 R_a/S_a、P/M 等标签，按原文命名系统保留。不能仅由中心手性标记假设轴手性已编码。
- `assignment`：绝对/相对构型、作者的测定/归属方法和不确定性；图中楔线经结构推导成 R/S 时写明推导，不能当成原文明确命名。
- `mixture`：外消旋、对映体/非对映体/几何异构体混合物及组成。未知组成不填 1:1；外消旋不是 achiral。
- `chirality_status`：specified / partially_specified / unspecified / racemic / mixture / achiral / unresolved / not_applicable。
- `geometry_status`：specified / partially_specified / unspecified / mixture / not_applicable / unresolved。
- `representation_status`：encoded / partial / not_encoded / not_applicable / unresolved。图上已指定而工具无法表达的轴/环立体信息保留在其他字段并用 partial/not_encoded，不能丢掉后说已完整编码。
- `isomeric_smiles`：与主 `smiles` 保持同一立体图的异构 SMILES；确定的四面体中心使用 @/@@，适用双键保留 / 和反斜杠。结构/构型未定时明确留空或保留已知部分，不能编造。
- `identity_key`：可选的工具生成 stereo_graph_key，仅作异构分子图比较，不是样品唯一身份。不足以区分的混合物、外消旋体与未知构型必须另建 molecule_id。

每个分子必须填写上述状态事实，无法判断就 unresolved 并给出证据/原因；不能全部用空骨架替代立体检查。achiral/not_applicable 需要结构依据，以 structure_derivation/inference 记录，不从缺少标记推断。分别命名的 E/Z 或 R/S 异构体分别建档；一个只给比例的混合物保留混合物身份，不能冒充单一结构或凭比例创建未证实的绝对构型。

## 工具与实验结果

`python scripts/reaction_tools.py standard stereo --smiles 'C/C=C/C'`

工具返回规范异构 SMILES、保留立体差异的图键、潜在/已指定中心和双键。输出只是结构计算，保存为派生事实时需原图证据与 derivation。RDKit 的轴手性/环立体支持有限；无输出不等于无手性。禁止用去立体的 connectivity_smiles 合并分子。

ee、er、dr、E/Z、cis/trans 比例分别作为 `outcomes[kind=selectivity]`，每项有单独 target_molecule_ids、value、scope（指标名、异构体顺序与分子 ID）、measurement_basis 及证据。保留 >、<、范围、主/次构型和检测方法；99% ee 不意味着 99:1 er，收率不是选择性。仅有 ee 不倒推出绝对构型；未分配的主异构体不猜 R/S/E/Z。选择性是该次实验的结果，不能无依据推广到同一产物的其他实验。

导入后复查：相反构型是否不同分子 ID，图上楔线及双键是否与 SMILES 相符，原文未定构型是否仍未定，各选择性结果的对象/分母/指标/顺序是否明确。所有机器结果仍 needs_review。
