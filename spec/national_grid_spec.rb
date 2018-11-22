require_relative '../spec/spec_helper'

RSpec.describe NationalGrid do
  context 'When given null values' do
    it 'throws an error' do
      expect(OSGB36toWGS84(531061, 104458)).to eq(50.824841, -0.14057994)
    end
  end
end
