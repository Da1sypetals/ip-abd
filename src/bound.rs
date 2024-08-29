pub struct Boundary;
impl Boundary {
    // Left Right; Bottom Top
    pub fn left_bot() -> glm::Vec2 {
        return glm::vec2(-1., -1.);
    }

    pub fn left_top() -> glm::Vec2 {
        return glm::vec2(-1., 1.);
    }

    pub fn right_bot() -> glm::Vec2 {
        return glm::vec2(1., -1.);
    }

    pub fn right_top() -> glm::Vec2 {
        return glm::vec2(1., 1.);
    }
    pub fn points() -> Vec<glm::Vec2> {
        vec![
            Boundary::left_bot(),
            Boundary::left_top(),
            Boundary::right_bot(),
            Boundary::right_top(),
        ]
    }

    pub fn bot() -> (glm::Vec2, glm::Vec2) {
        return (Self::left_bot(), Self::right_bot());
    }

    pub fn top() -> (glm::Vec2, glm::Vec2) {
        return (Self::left_top(), Self::right_top());
    }

    pub fn left() -> (glm::Vec2, glm::Vec2) {
        return (Self::left_bot(), Self::left_top());
    }

    pub fn right() -> (glm::Vec2, glm::Vec2) {
        return (Self::right_bot(), Self::right_top());
    }

    // to avoid some numerical issues
    pub fn edges() -> Vec<(glm::Vec2, glm::Vec2)> {
        vec![
            (
                glm::vec2(Boundary::left_bot().x * 1.1, Boundary::left_bot().y),
                glm::vec2(Boundary::right_bot().x * 1.1, Boundary::right_bot().y),
            ),
            (
                glm::vec2(Boundary::left_top().x * 1.1, Boundary::left_top().y),
                glm::vec2(Boundary::right_top().x * 1.1, Boundary::right_top().y),
            ),
            (
                glm::vec2(Boundary::left_bot().x, Boundary::left_bot().y * 1.1),
                glm::vec2(Boundary::left_top().x, Boundary::left_top().y * 1.1),
            ),
            (
                glm::vec2(Boundary::right_bot().x, Boundary::right_bot().y * 1.1),
                glm::vec2(Boundary::right_top().x, Boundary::right_top().y * 1.1),
            ),
        ]
    }
}
