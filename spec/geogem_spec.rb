require_relative 'spec_helper'

RSpec.describe Point do
  describe ' #distance' do
    it 'returns distance from two points' do
      expect(Point.new("400000 200000").distance(Point.new("400100 200000")))
      .to eq(100)
    end
  end

  describe ' #point_buff' do
    it 'returns circle geometry buffering point' do
      expect(Point.new("500000 150000").point_buff(100,12))
      .to eq("500100.0 150000.0, 500086.60254037846 150050.0, 500050.0 150086.60254037843, 500000.0 150100.0, 499950.0 150086.60254037843, 499913.39745962154 150050.0, 499900.0 150000.0, 499913.39745962154 149950.0, 499950.0 149913.39745962157, 500000.0 149900.0, 500050.0 149913.39745962157, 500086.60254037846 149950.0, 500100.0 150000.0")
    end
  end

  describe ' #to_wkt' do 
    it 'returns the point in Well Known Text format' do 
      expect(Point.new("200000 300000").to_wkt)
      .to eq("POINT(200000 300000)")
    end
  end
end

RSpec.describe Polygon do 
  describe ' #centre_point' do
    it 'returns centre point of a polygon' do
      expect(Polygon.new("520575 170388, 520617 170405, 520624 170389, 520582 170371, 520575 170388").centre_point)
      .to eq("520599.5 170388.0")
    end
  end

  describe ' #bbox' do
    it 'returns bounding box geometry of a polygon' do
      expect(Polygon.new("520575 170388, 520617 170405, 520624 170389, 520582 170371, 520575 170388").bbox)
      .to eq("520575 170371, 520624 170405")
    end
  end

  describe ' #area' do
    it 'returns the area of a polyon' do
      expect(Polygon.new("0 0, 10 0, 10 10, 0 10, 0 0").area)
      .to eq(100)
    end
  end

  describe ' #intersect?' do
    it 'returns true if 2 line segments intersect' do
      my_obj = Polygon.new("0 0")
      expect(my_obj.send(:intersect?, [Point.new("0 0"), Point.new("2 2")], [Point.new("2 0"), Point.new("0 2")]))
      .to eq (true)
    end
  end

  describe ' #self_intersects?' do
    it 'returns true if a polygon intersects with it\'s self' do 
      expect(Polygon.new("0 0, 2 2, 0 2, 2 0, 0 0").self_intersects?)
      .to eq(true)
    end
  end

  describe '#to_wkt(polygon)' do
    it 'returns the wkt of a polygon' do
      expect(Polygon.new("Polygon((0 0, 10 0, 10 10, 0 10, 0 0))").to_wkt)
      .to eq("POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))")
    end
  end

  describe '#intersects?' do
    it 'returns true if 2 polygons intersect' do 
      expect(Polygon.new("0 0, 10 0, 10 10, 0 10, 0 0").intersects?(Polygon.new("5 5, 20 0, 20 30, 0 20, 5 5")))
      .to eq(true)
    end
  end

  describe '#point_in_poly?' do 
    it 'returns true if a point lies inside a polygon' do
      expect(Polygon.new("0 0, 100 0, 100 100, 0 100, 0 0").send(:point_in_poly?, Point.new("50 50")))
      .to eq(true)
    end
  end
end