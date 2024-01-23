#include "imagingtable.h"

#include <aocommon/logger.h>
#include <schaapcommon/facets/facet.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <map>
#include <memory>

using aocommon::Logger;

ImagingTable::ImagingTable(const std::vector<EntryPtr>& entries)
    : _entries(entries), _independentGroups(), _facets(), _squaredGroups() {
  Update();
}

ImagingTable::ImagingTable(
    const ImagingTable& other,
    std::function<bool(const ImagingTableEntry&)> isSelected)
    : _entries(), _independentGroups(), _facets(), _squaredGroups() {
  std::copy_if(other._entries.begin(), other._entries.end(),
               std::back_inserter(_entries),
               [&](const EntryPtr& entry) { return isSelected(*entry); });
  Update();
}

void ImagingTable::Print() const {
  Logger::Info << "=== IMAGING TABLE ===\n"
                  "       # Pol Ch JG ²G FG FI In Freq(MHz)\n";
  for (size_t i = 0; i != _independentGroups.size(); ++i) {
    Logger::Info << "| Independent group:\n";
    const ImagingTable independent(_independentGroups[i]);

    const ImagingTable::Groups& squaredGroups = independent.SquaredGroups();
    for (size_t s = 0; s != squaredGroups.size(); ++s) {
      const ImagingTable::Group& squared = squaredGroups[s];

      for (size_t e = 0; e != squared.size(); ++e) {
        if (s == 0 && e == 0)
          Logger::Info << "+-";
        else if ((i + 1) == _independentGroups.size())
          Logger::Info << "  ";
        else
          Logger::Info << "| ";

        if (e == 0)
          Logger::Info << "+-";
        else if ((s + 1) == squaredGroups.size())
          Logger::Info << "  ";
        else
          Logger::Info << "| ";

        PrintEntry(*squared[e]);
      }
    }
  }
}

void ImagingTable::PrintEntry(const ImagingTableEntry& entry) {
  std::ostringstream str;

  str << std::setw(2) << entry.index << "  ";
  str << aocommon::Polarization::TypeToShortString(entry.polarization) << "  ";
  str << std::setw(2) << entry.outputChannelIndex << " ";
  str << std::setw(2) << entry.joinedGroupIndex << " ";
  str << std::setw(2) << entry.squaredDeconvolutionIndex << " ";
  str << std::setw(2) << entry.facetGroupIndex << " ";
  str << std::setw(2) << entry.facetIndex << " ";
  str << std::setw(2) << entry.outputIntervalIndex << "  ";
  str << round(entry.bandStartFrequency * 1e-6) << "-"
      << round(entry.bandEndFrequency * 1e-6) << " (" << entry.inputChannelCount
      << ")";

  Logger::Info << "J-" << str.str() << '\n';
}

void ImagingTable::updateGroups(
    Groups& groups, std::function<size_t(const ImagingTableEntry&)> getIndex,
    std::function<bool(const ImagingTableEntry&)> isSelected) const {
  std::map<size_t, Group> groupMap;

  for (const EntryPtr& e : _entries) {
    if (isSelected(*e)) {
      groupMap[getIndex(*e)].push_back(e);
    }
  }

  groups.clear();
  for (auto& item : groupMap) {
    groups.emplace_back(std::move(item.second));
  }
}

void ImagingTable::AssignGridDataFromPolarization(
    aocommon::PolarizationEnum polarization) {
  for (Group& group : _squaredGroups) {
    const EntryPtr& sourceEntry = *std::find_if(
        group.begin(), group.end(), [polarization](const EntryPtr& e) {
          return e->polarization == polarization;
        });
    for (EntryPtr& entryPtr : group) {
      if (entryPtr != sourceEntry) {
        entryPtr->AssignGridData(*sourceEntry);
      }
    }
  }
}

std::unique_ptr<radler::WorkTable> ImagingTable::CreateDeconvolutionTable(
    int n_deconvolution_channels, CachedImageSet& psf_images,
    CachedImageSet& model_images, CachedImageSet& residual_images) const {
  if (_entries.empty()) {
    throw std::runtime_error(
        "Can not create a DeconvolutionTable from an empty ImagingTable.");
  }

  if (_facets.empty()) {
    throw std::runtime_error(
        "Member _facets of ImagingTable is empty. Probably Update() was not "
        "called after adding entries.");
  }

  // In a DeconvolutionTable the output channel indices range from
  // 0 to (#channels - 1). In an ImagingTable that forms an indepent group,
  // output channel indices may start at a higher index.

  // Assume that the first entry has the lowest index and that the last entry
  // has the highest output channel index.
  const int channel_index_offset = _entries.front()->outputChannelIndex;
  const int n_original_channels =
      _entries.back()->outputChannelIndex + 1 - channel_index_offset;

  // Gather dd-psf info from ImagingTable
  // The ImagingTable is a fairly free format
  // The following addional assumptions are needed to create a valid
  // WorkTable with ddp-psfs
  //
  // * The dd-psfs are the same, and in the same order for all images
  // * the facet_ids of the dd-psfs form a simple sequence 0,1,2,...
  //

  std::vector<radler::PsfOffset> psf_offsets;
  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> psf_facets;

  // There could be multiple independent groups so the first
  // facetGroupIndex is not necessarily zero.
  // Note that looping over _facetGroups.front() is not possible here, since it
  // does not contain the DD PSF entries.
  size_t first_facet_group_index = (*_entries.begin())->facetGroupIndex;

  // Although an ImagingTable may contain entries for multiple polarizations,
  // there should only be DD PSF entries for the first polarization, since the
  // DD PSF layout is equal for all polarizations.
  for (const EntryPtr& entry_ptr : _entries) {
    if (entry_ptr->facetGroupIndex != first_facet_group_index) break;
    if (!entry_ptr->isDdPsf) continue;
    const schaapcommon::facets::PixelPosition centre =
        entry_ptr->facet->GetTrimmedBoundingBox().Centre();
    psf_offsets.emplace_back(centre.x, centre.y);
    psf_facets.push_back(entry_ptr->facet);
  }

  auto table = std::make_unique<radler::WorkTable>(
      std::move(psf_offsets), n_original_channels, n_deconvolution_channels,
      channel_index_offset);
  int max_squared_index = -1;

  for (const EntryPtr& entry_ptr : _facets.front()) {
    assert(entry_ptr);

    if (entry_ptr->imageCount >= 1) {
      CachedImageSet* psf_images_ptr = nullptr;

      // Only set psf_images_ptr for the first entry of each squared group.
      // This way, CreateDeconvolutionEntry() only creates psf accessors for
      // those entries.
      if (int(entry_ptr->squaredDeconvolutionIndex) > max_squared_index) {
        max_squared_index = entry_ptr->squaredDeconvolutionIndex;
        psf_images_ptr = &psf_images;
      }

      std::unique_ptr<radler::WorkTableEntry> real_entry =
          entry_ptr->CreateDeconvolutionEntry(
              channel_index_offset, psf_images_ptr, model_images,
              residual_images, psf_facets, false);
      table->AddEntry(std::move(real_entry));
    }

    if (entry_ptr->imageCount == 2) {
      std::unique_ptr<radler::WorkTableEntry> imaginary_entry =
          entry_ptr->CreateDeconvolutionEntry(channel_index_offset, nullptr,
                                              model_images, residual_images,
                                              psf_facets, true);
      table->AddEntry(std::move(imaginary_entry));
    }
  }
  return table;
}
