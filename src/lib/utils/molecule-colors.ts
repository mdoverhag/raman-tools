// Consistent color mapping for molecules
const moleculeColors: Record<string, string> = {
  // Raman molecules - all blue
  DTNB: "bg-blue-400/10 text-blue-400",
  MBA: "bg-blue-400/10 text-blue-400",
  TFMBA: "bg-blue-400/10 text-blue-400",

  // Target molecules - all pink
  IgG: "bg-pink-400/10 text-pink-400",
  BSA: "bg-pink-400/10 text-pink-400",
  HER2: "bg-pink-400/10 text-pink-400",
  EpCAM: "bg-pink-400/10 text-pink-400",
  TROP2: "bg-pink-400/10 text-pink-400",

  // Default for unknown molecules
  default: "bg-gray-400/10 text-gray-400",
};

export function getMoleculeClasses(molecule: string): string {
  return moleculeColors[molecule] || moleculeColors.default;
}
