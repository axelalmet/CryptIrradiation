#include "StromalCellProliferativeType.hpp"

StromalCellProliferativeType::StromalCellProliferativeType()
    : AbstractCellProliferativeType(4)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StromalCellProliferativeType)
