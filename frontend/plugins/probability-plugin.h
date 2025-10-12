#pragma once

#include "../imgui-renderer.h"

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace math_gui::plugins
{

enum class DistributionType
{
	Continuous,
	Discrete
};

using PdfFunction = std::function<float(float)>;
using PdfFactory = std::function<PdfFunction(const std::vector<float>&)>;
using CdfFunction = std::function<float(float)>;
using CdfFactory = std::function<CdfFunction(const std::vector<float>&)>;
using PpfFunction = std::function<float(float)>;
using PpfFactory = std::function<PpfFunction(const std::vector<float>&)>;
using StatisticValue = std::optional<float>;
using StatisticFactory = std::function<StatisticValue(const std::vector<float>&)>;

struct StatisticDefinition
{
	std::string label;
	std::string formula;
	StatisticFactory compute;
};

struct DistributionDefinition
{
	std::array<float, 2> domain;
	std::array<float, 2> renderDomain;
	DistributionType type;
	std::vector<float> parameters;
	std::vector<std::array<float, 2>> parameterDomains;
	std::vector<std::string> parameterNames;
	std::vector<std::string> parameterRanges;
	std::vector<std::string> parameterDescriptions;
	PdfFactory makePdf;
	std::vector<bool> parameterIsIntegral;
	CdfFactory makeCdf;
	PpfFactory makePpf;
	std::vector<StatisticDefinition> statistics;
};

struct DistributionEntry
{
	std::string label;
	DistributionDefinition definition;
};

void RenderProbabilityWindow(ImGuiRenderer::FrameState& state);

} // namespace math_gui::plugins
