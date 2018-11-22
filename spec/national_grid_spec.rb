require_relative '../spec/spec_helper'

RSpec.describe NationalGrid do
  context 'When given null values' do
    let(:easting) { }
    let(:northing) { }

    it 'throws an error' do
      xyz
    end
  end

  context 'When given Easting Northing' do
    let(:easting) { }
    let(:northing) { }

    it 'returns the correct Lat Long' do
      xyz
    end
  end
end
